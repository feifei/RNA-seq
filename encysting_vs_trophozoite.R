library(reshape2)
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtools)
library(edgeR)


get_logFC.sig <- function(name = "human", ...) {
  ##############
  # Make a table of logFC with genes has logFC > 1 in any comparisons
  ##############
  logFC_trophozoite_infile <- paste("data/", name, ".", "trophozoite.logFC.rda", sep = "")
  load(logFC_trophozoite_infile)
  trophozoite_logFC_mat <- logFC_mat

  logFC_encysting_infile <- paste("data/", name, ".", "encysting.logFC.rda", sep = "")
  load(logFC_encysting_infile)
  encysting_logFC_mat <- logFC_mat
  

  colnames(trophozoite_logFC_mat) <- paste("trophozoite", colnames(trophozoite_logFC_mat), sep = "_")
  colnames(encysting_logFC_mat) <- paste("encysting", colnames(encysting_logFC_mat), sep = "_")
  
  # merge by row.names, remove first column Row.names
  logFC_mat <- merge(trophozoite_logFC_mat, encysting_logFC_mat, by=0)
  logFC_mat.sig <- logFC_mat[apply(logFC_mat[,-1], MARGIN = 1, function(x) any(x >= 1)), ]
  rownames(logFC_mat.sig) <- logFC_mat.sig$Row.names
  logFC_mat.sig <- logFC_mat.sig[,-1]
  
  if (name == "human") {
    load('data/genes_go_kegg.rda') # genes
    genes_sub <- genes[genes$ensembl_gene_id %in% rownames(logFC_mat.sig), ]
    description <- genes_sub$description[match(rownames(logFC_mat.sig), genes_sub$ensembl_gene_id)]
  } else {
    annotation_tab <- "data/wb_annotation.tab"
    annotation <- read.csv(file=annotation_tab, head=TRUE, sep="\t")
    description <- annotation$Description[match(rownames(logFC_mat.sig), annotation$Geneid)]
  }
  
  logFC_mat.sig$description <- description
  # Put description in second column
  logFC_mat.sig <- logFC_mat.sig[,c(ncol(logFC_mat.sig), 1:ncol(logFC_mat.sig)-1)]
  logFC_mat.sig <- logFC_mat.sig[order(logFC_mat.sig$trophozoite_Infection_1.5vs0, decreasing = TRUE), ]

  outfile <- paste("data/", name, ".", "encysting_trophozoite.logFC.csv", sep = "")
  write.csv(logFC_mat.sig, file = outfile)
  return(logFC_mat.sig)
}


get_logCPM.sig <- function(name = "human", logFC_mat.sig, ...) {
  #################
  # Get the logCPM, normalized with all the data combined
  # Similar to the ones in DE analysis
  ###################
  meta_file <- "data/metadata2.tab"
  d <- read.csv(meta_file, sep = ",")
  count_files <- paste(d$identifier, ".", name, ".counts", sep = "")
  d$count_files <- count_files
  if (name == "human") {
    d <- d[-c(1, 6, 11, 17, 21), ] # Remove Giairdia controls
  } else {
    d <- d[-c(2, 7, 12, 16, 25), ] # Remove Human controls
  }
  d <- d[with(d, order(hour, condition)), ]
  
  group <- factor(paste(d$cell, d$condition, d$hour, sep = "_"))
  d <- cbind(d, group = group)
  y <- readDGE(d$count_files, path = "counts", group = d$group, header = FALSE)
  colnames(y$counts) <- group
 
  if (name == "human") {
    # rownames, aka geneid has suffix which were trimmed, and duplciated entry with _PAR_Y were removed
    old_rownames <- rownames(y$counts)
    bad_row_idx <- grep("_PAR_Y", old_rownames)
    sub_rownames <- old_rownames[- bad_row_idx]
    new_rownames <- sapply(strsplit(sub_rownames, split = '.', fixed = TRUE), function(x) (x[1]))
    y$counts <- y$counts[- bad_row_idx, ]
    rownames(y$counts) <- new_rownames
    
    load('data/genes_go_kegg.rda') # genes
    y$genes <- genes[match(rownames(y$counts), genes$ensembl_gene_id), ]
  } else {
    annotation_tab <- "data/wb_annotation.tab"
    annotation <- read.csv(file=annotation_tab, head=TRUE, sep="\t")
    sub_annotation <- annotation[match(rownames(y$counts), annotation$Geneid), ]
    y$genes <- sub_annotation
  }

  CPM <- cpm(y)
  keep <- rowSums(CPM > 2) >= 3
  table(keep)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  logCPM <- cpm(y, log = TRUE, prior.count = 2)

  ###########
  # Make a table of logCPM with genes that is DE in any comparison
  ###########
  logCPM.sig <- logCPM[rownames(logCPM) %in% rownames(logFC_mat.sig),]
  logCPM.sig <- as.data.frame(logCPM.sig)
  logCPM.sig <- logCPM.sig[match(rownames(logFC_mat.sig), rownames(logCPM.sig)), ] # Order by logFC_mat
  logCPM.sig$description <- logFC_mat.sig$description
  # Reording changing the colnames to unique ones, which makes the hours wrong for the plot
  logCPM.sig <- logCPM.sig[,c(ncol(logCPM.sig), 1:ncol(logCPM.sig)-1)]
  return(logCPM.sig)
}



grid_arrange_shared_legend <- function(plots) {
  # plots is a list of ggplots
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

plot_logCPM.sig <- function(name = "human", logCPM.sig, ...) {
  pdf(file = paste("plots/", name, ".encysting_trophozoite.exp_profiles.pdf", sep=""), height = 10, width = 7.5)
  
  p <- list()
  i <- 0
  for (geneid in rownames(logCPM.sig)) {
    i <- i + 1
    row <- logCPM.sig[i, ]
    row_m <- melt(row, value.name="logCPM")
    row_m[c("cell", "condition", "hour")] <- str_split_fixed(row_m$variable, "_", 3)
    row_m$hour <- gsub("\\.\\d$", "", row_m$hour) # need to fix 1 and 4 
    row_m$hour[row_m$hour == 4] <- 4.5
    row_m$hour[row_m$hour == 1] <- 1.5
    row_m$hour <- factor(row_m$hour, levels=c(0, 1.5, 3, 4.5))
    row_m$cell <- factor(row_m$cell, levels=c("trophozoite", "encysting"))
    description = row$description
    p[[i]] <- ggplot(row_m, aes(x=hour, y=logCPM, color=cell, group=cell)) + geom_point() +
      stat_smooth(se=FALSE,method="loess") + labs(title = paste(geneid, description, sep="\n")) +
      theme(plot.title = element_text(size=7)) + scale_color_manual(values=c("#999999", "#F8766D"))
    if (i == 6) {
      grid_arrange_shared_legend(p)
      i <- 0
      p <- list()
    }
  }
  
  dev.off()
}

logFC.sig.human <- get_logFC.sig("human")
logCPM.sig.human <- get_logCPM.sig("human", logFC_mat.sig = logFC.sig.human)
write.csv(logFC.sig.human, file = "data/human.encysting_trophozoite.logFC.csv")
write.csv(logCPM.sig.human, file = "data/human.encysting_trophozoite.logCPM.csv")
plot_logCPM.sig("human", logCPM.sig.human)

logFC.sig.giardia <- get_logFC.sig("giardia")
logCPM.sig.giardia <- get_logCPM.sig("giardia", logFC_mat.sig = logFC.sig.giardia)
write.csv(logFC.sig.giardia, file = "data/giardia.encysting_trophozoite.logFC.csv")
write.csv(logCPM.sig.giardia, file = "data/giardia.encysting_trophozoite.logCPM.csv")
plot_logCPM.sig("giardia", logCPM.sig.giardia)
