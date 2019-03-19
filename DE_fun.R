
library(edgeR)
library(DT)
library(GO.db)
library(goseq)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(org.Hs.eg.db)
library(limma)


qlf_dictionary <- function(fit, my.contrasts, ...) {
  # To pair up contrast and qlf
  qlf_d <- list()
  for (contrast in dimnames(my.contrasts)$Contrasts) {
    qlf <- glmQLFTest(fit, contrast=my.contrasts[, contrast])
    qlf_d <- c(qlf_d, list(name = qlf))
  }
  names(qlf_d) <- dimnames(my.contrasts)$Contrasts
  return(qlf_d)
}

plot_logFCvslogCPM <- function(qlf, contrast, ...) { 
  cat("\n### Log-fold change against log-counts per million\n")
  cat("The test results are visualized in te following smear plot. 
      Genes that are significantly DE with an FDR of 5% are highlighted in red and blue.\n\n")
  plotMD(qlf, main=contrast)
  abline(h=c(-1,1), col="blue")
}

plot_pheatmap <- function(qlf, summ, contrast, tag = "", ...) {
  cat("\n\n### Heatmap of the moderated log-counts-per-million\n")
  cat("\nHeatmap of the moderated log-counts-per-million of the DE genes in this comparison.")
  size <- sum(summ[c(1,3)])
  if (size == 0) return()
  DE.top <- topTags(qlf, n=size, sort.by="logFC")
  top.genes <- row.names(DE.top)
  top.logcpm <- logCPM[top.genes, ]
  pheatmap(top.logcpm, cluster_cols = FALSE,
           show_rownames=FALSE, color=diverge_hcl(100), main=paste(contrast, "log-counts-per-million", tag))
}

plot_valcano <- function(all.genes.t, contrast, ...) {
  cat("\n\n### Volcano plot\n")
  cat("A volcano plot the fold change on the x-axis and the statistial significance on the y-axis
      (the -log10 of the p-value. Genes that are highly dysregulated are farther to the left and right sides, 
      while highly significant changes appear higher on the plot. 
      Genes with FDR<0.05 are colored red, and abs(logFC)>1 are colored orange, and green if both.\n\n")
  with(all.genes.t, plot(logFC, -log10(PValue), pch=20, main=paste("Volcano plot", contrast)))
  with(subset(all.genes.t, FDR<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
  with(subset(all.genes.t, abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="orange"))
  with(subset(all.genes.t, FDR<.05 & abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="green"))
}

export_DE <- function(DE.genes.t, directory, name, cell, contrast, tag = "", ...) {
  cat("\n\n### Export significant cases\n")
  cat("The table lists all the DE with abs(logFC) > 1 and FDR< 0.05. \n\n")
  
  colSize = ncol(DE.genes.t)
  # Workaround in printing table
  print(
    htmltools::tagList(
      datatable(DE.genes.t[abs(DE.genes.t$logFC) > 1, ], rownames = FALSE, 
                options = list(scrollX = TRUE)) %>%
        formatStyle(columns = c(1:colSize), fontSize = '80%') %>%
        formatSignif(c(8:colSize), 2) %>%
        formatRound(c(5:7), 2)
    )
  )
  
  DE.genes.t <- format(DE.genes.t, trim=TRUE, digits=2, nsmall=2)
  write.csv(DE.genes.t, file=paste(directory, name, ".", cell, ".", contrast, ".sig_results", tag, ".csv", sep=""), 
            row.names = FALSE)
}

export_all <- function(all.genes.t, directory, name, cell, contrast, tag = "", ...) {
  cat("\n\n### Export all genes\n")
  cat("The table lists all the genes. \n\n")
  
  all.genes.t <- format(all.genes.t, trim=TRUE, digits=2, nsmall=2)
  write.csv(all.genes.t, file=paste(directory, name, ".", cell, ".", contrast, ".all", tag, ".csv", sep=""), 
            row.names = FALSE)
}

export_cpm <- function(DE.genes.t, directory, name, cell, contrast, tag = "", ...) {
  cat("\n\n### Export logCPM values of the DE\n")
  cat("Inspect the depth-adjusted reads per million for the top differentially expressed. \n\n")
  
  DE.cpm <- logCPM[rownames(DE.genes.t),]
  if (nrow(DE.genes.t) == 1) {
    # numeric instead of matrix
    DE.cpm <- matrix(DE.cpm, nrow=1)
    colnames(DE.cpm) <- group
  }  
  
  colSize <- ncol(DE.cpm)
  print(
    htmltools::tagList(
      datatable(DE.cpm, options = list(scrollX = TRUE)) %>%
        formatStyle(columns = c(1:colSize), fontSize = '80%') %>%
        formatRound(c(1:colSize), 1)
    )
  ) 
  
  DE.cpm <- format(DE.cpm, digits=2, nsmall=1)
  write.csv(DE.cpm, file=paste(directory, name, ".", cell, ".", contrast, ".sig_results.cpm", tag, ".csv", sep=""))
}

go_enrich_plot <- function(GOmap, genes, lengthData, contrast, tag = "", extra_tag = "", ...) {
  # Find GO enrichment, and plot term vs num
  pwf=nullp(genes, bias.data=lengthData, plot.fit=FALSE)
  GO.wall=goseq(pwf,gene2cat=GOmap)
  
  GO.wall.sub = GO.wall[GO.wall$numDEInCat >= 5 | p.adjust(GO.wall$over_represented_pvalue, method="BH")<1,]
  if (nrow(GO.wall.sub) <= 1) return()
  GO.wall.sub$term <- reorder(GO.wall.sub$term, GO.wall.sub[-1]$numDEInCat)
  GO.wall.sub.m <- melt(GO.wall.sub, id=c("category","term", "over_represented_pvalue","under_represented_pvalue","ontology"))
  
  # GO.wall.sub$ontology[is.na(GO.wall.sub$ontology)] <- "NA"
  GO.wall.sub.sub <- head(GO.wall.sub[order(GO.wall.sub$numDEInCat,decreasing=TRUE),], n = 50)
  p <- ggplot(subset(GO.wall.sub.sub, !is.na(ontology)), 
              aes(term, numDEInCat, fill=ontology))
  # p <- ggplot(GO.wall.sub, aes(term, numDEInCat, fill=ontology))
  p <- p + geom_bar(stat="identity", position=position_dodge()) +
    coord_flip() + labs(title=paste(contrast, extra_tag, tag)) +
    scale_x_discrete(label=function(x) substr(x, 1, 40))
  return (list("GO.wall" = GO.wall, "GO.wall.sub.m" = GO.wall.sub.m, "plot" = p))
}

go_up_down <- function(GOmap, up_down, genes, lengthData, contrast, extra_tag = "", tag = "", ...) {
  sub_up_down <- up_down[up_down %in% names(lengthData)]
  sub_genes <- as.integer(names(genes) %in% sub_up_down)
  names(sub_genes) <- names(genes)
  g <- go_enrich_plot(GOmap, sub_genes, lengthData, contrast, extra_tag = extra_tag, tag = tag)
  if (is.null(g)) return()
  p <- g$plot
  GO.wall = g$GO.wall
  GO.wall.sub.m = g$GO.wall.sub.m
  if (length(unique(GO.wall.sub.m$ontology)) == 3) {
    p <- p + scale_fill_manual(values=c("#F8766D", "#00BA38", "#999999"))
  } else {
    p <- p + scale_fill_manual(values=c("#F8766D", "#00BA38", "#619CFF", "#999999"))
  }
  print (p)
  cat('\r\n\r\n')
  
  GO.wall$FDR <- p.adjust(GO.wall$over_represented_pvalue, method="BH")
  GO.wall.sig <- GO.wall[GO.wall$FDR < 0.05, ]
  if (nrow(GO.wall.sig) == 0) return()
  write.csv(GO.wall.sig, file=paste(directory, name, ".", cell, ".", contrast,
                                    ".sig_results.GO.", extra_tag, tag, ".csv", sep=""),
            row.names = FALSE)
  
  enriched.GO=GO.wall.sig$category
  cat("------------ Enriched", extra_tag, "-regulated GO ------------\r\n")
  
  # Workaround in printing table in loop
  print(
    htmltools::tagList(
      datatable(GO.wall.sig, options = list(scrollX = TRUE), rownames = FALSE) %>%
        formatStyle(columns = c(1:colSize), fontSize = '80%')
    )
  )
  
  print(enriched.GO)
  cat("\n")
  for(go in enriched.GO){
    goterm <- GOTERM[[go]]
    if (!is.null(goterm)) {
      cat("GOID: ", GOID(goterm), "\n")
      cat("Term: ", Term(goterm), "\n")
      cat("Ontology: ", Ontology(goterm), "\n")
      cat("Definition: ", Definition(goterm), "\n")
      cat("--------------------------------------\n\n") 
    }
  }
}

go_analysis <- function(y, GOmap, isDE, dt, contrast, tag = "", ...) {
  ########################
  # Go enrichment analysis with goseq, which adjust for gene length bias
  ########################
  cat("\n\n\n\n### GO analysis with goseq \n\n")
  isDE_int <- as.integer(isDE)
  names(isDE_int) <- rownames(y)
  
  isDE_int <- isDE_int[!is.na(isDE_int)]
  
  length_vect <- y$genes$length
  names(length_vect) <- rownames(y)
  g <- go_enrich_plot(GOmap, isDE_int, length_vect, contrast)
  if (is.null(g)) return()
  p <- g$plot
  print (p)
  GO.wall.sub.m = g$GO.wall.sub.m
  
  # Plot ontology big class 
  p <- ggplot(GO.wall.sub.m, aes(ontology, value, fill=variable))
  print(p + geom_bar(stat="identity", position=position_dodge()) +
          coord_flip() + labs(title=paste(contrast, tag)))
  
  
  ####
  # Separate Up and Down regulation
  up = rownames(y)[which(dt == 1)]
  down = rownames(y)[which(dt == -1)]
  if (length(up) > 1){
    go_up_down(GOmap, up, isDE_int, length_vect, contrast, extra_tag = "UP")
  }
  if (length(down) > 1) {
    go_up_down(GOmap, down, isDE_int, length_vect, contrast, extra_tag = "DOWN")
  }
}


kegg_up_down <- function(up_down, lengthData, isDE_int, KEGGmap, extra_tag = "UP", ...) {
  up_down <- up_down[up_down %in% names(lengthData)]
  isDE_int <- isDE_int[names(isDE_int) %in% up_down]

  pwf=nullp(isDE_int, bias.data=lengthData, plot.fit=FALSE)
  KEGG=goseq(pwf,gene2cat=KEGGmap)
  
  KEGG$FDR <- p.adjust(KEGG$over_represented_pvalue, method="BH")
  KEGG_sig <- KEGG[KEGG$FDR < 0.05, ]
  if (nrow(KEGG_sig) == 0) return()
  write.csv(KEGG_sig, file=paste(directory, name, ".", cell, ".", contrast,
                                    ".sig_results.KEGG.", extra_tag, tag, ".csv", sep=""),
            row.names = FALSE)
  
  # Workaround in printing table in loop
  print(
    htmltools::tagList(
      datatable(KEGG_sig, options = list(scrollX = TRUE), rownames = FALSE) %>%
        formatStyle(columns = c(1:colSize), fontSize = '80%')
    )
  )
}

kegg_analysis <- function(up_down, lengthData, isDE_int, KEGGmap, ...) {
  ########################
  # KEGG enrichment analysis with goseq, which adjust for gene length bias
  ########################
  cat("\n\n\n\n### KEGG analysis with goseq \n\n")
  isDE_int <- as.integer(isDE)
  names(isDE_int) <- rownames(y)
  isDE_int <- isDE_int[!is.na(isDE_int)]
  
  length_vect <- y$genes$length
  names(length_vect) <- rownames(y)
  
  up = rownames(y)[which(dt == 1)]
  down = rownames(y)[which(dt == -1)]
  
  kegg_up_down(up, lengthData, isDE_int, KEGGmap, extra_tag = "UP")
  kegg_up_down(down, lengthData, isDE_int, KEGGmap, extra_tag = "DOWN")
  
  }



########## Results from the following method is not good
go_kegg_up_down <- function(up_down, go_kegg_tag = "GO", extra_tag = "UP", ...) {
  if (go_kegg_tag == "GO") {
    go_kegg <- goana(up_down, species = "Hs")
    top_go_kegg <- topGO(go_kegg)
  } else if (go_kegg_tag == "KEGG") {
    go_kegg <- kegga(up_down, species = "Hs")
    top_go_kegg <- topKEGG(go_kegg)
  }
  go_kegg <- go_kegg[order(go_kegg$P.DE), ] # Order
  go_kegg <- go_kegg[go_kegg$P.DE < 0.05, ] # Filter
  top_go_kegg$P.DE <- format(top_go_kegg$P.DE, nsmall=1, digits = 3)
  colSize = ncol(top_go_kegg)
  # Workaround in printing table in loop
  print(
    htmltools::tagList(
      datatable(top_go_kegg, options = list(scrollX = TRUE)) %>%
        formatStyle(columns = c(1:colSize), fontSize = '80%')
    )
  )
  
  if (nrow(go_kegg) > 1) { 
    write.csv(go_kegg, 
              file=paste(directory, name, ".", cell, ".", contrast, ".sig_results.", 
                         go_kegg_tag, ".", extra_tag, tag, ".csv", sep=""))
  }
}

go_kegg_analysis <- function(y, isDE, dt, directory, name, cell, contrast, tag="", ...) {
  ########################
  # Go and KEGG enrichment analysis with limma
  ########################
  cat("\n\n\n\n### GO & KEGG analysis \n\n")
  
  ####
  # Separate Up and Down regulation
  up = y$genes$entrezgene[which(dt == 1)]
  down = y$genes$entrezgene[which(dt == -1)]
  
  go_kegg_up_down(up, go_kegg_tag = "GO", extra_tag = "UP")
  go_kegg_up_down(down, go_kegg_tag = "GO", extra_tag = "DOWN")
  go_kegg_up_down(up, go_kegg_tag = "KEGG", extra_tag = "UP")
  go_kegg_up_down(down, go_kegg_tag = "KEGG", extra_tag = "DOWN")
}
