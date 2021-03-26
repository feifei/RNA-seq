
library(edgeR)
library(DT)
library(GO.db)
library(goseq)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(KEGGREST)


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
  if (nrow(DE.genes.t) >= 10){
    cat(paste(rep("\n", 11), collapse = ""))
  } else {
    cat(paste(rep("\n", nrow(DE.genes.t)), collapse = ""))
  }
  
  DE.genes.t <- format(DE.genes.t, trim=TRUE, digits=2, nsmall=2)
  write.csv(DE.genes.t, file=paste(directory, name, ".", cell, ".", contrast, ".sig_results", tag, ".csv", sep=""), 
            row.names = FALSE)
}

export_all <- function(all.genes.t, directory, name, cell, contrast, tag = "", ...) {
  cat("\n\n\n\n\n\n### Export all genes\n")
  cat("The table lists all the genes. \n\n")
  
  all.genes.t <- format(all.genes.t, trim=TRUE, digits=2, nsmall=2)
  write.csv(all.genes.t, file=paste(directory, name, ".", cell, ".", contrast, ".all", tag, ".csv", sep=""), 
            row.names = FALSE)
}

export_logcpm <- function(DE.genes.t, directory, name, cell, contrast, tag = "", ...) {
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
  if (nrow(DE.cpm) >= 10){
    cat(paste(rep("\n", 11), collapse = ""))
  } else {
    cat(paste(rep("\n", nrow(DE.cpm)), collapse = ""))
  }
  
  DE.cpm <- format(DE.cpm, digits=2, nsmall=1)
  write.csv(DE.cpm, file=paste(directory, name, ".", cell, ".", contrast, ".sig_results.logcpm", tag, ".csv", sep=""))
}

export_cpm <- function(DE.genes.t, directory, name, cell, contrast, tag = "", ...) {
  cat("\n\n### Export logCPM values of the DE\n")
  cat("Inspect the depth-adjusted reads per million for the top differentially expressed. \n\n")
  
  DE.cpm <- CPM[rownames(DE.genes.t),]
  if (nrow(DE.genes.t) == 1) {
    # numeric instead of matrix
    DE.cpm <- matrix(DE.cpm, nrow=1)
    colnames(DE.cpm) <- group
  }  
  
  DE.cpm <- format(DE.cpm, digits=2, nsmall=1)
  write.csv(DE.cpm, file=paste(directory, name, ".", cell, ".", contrast, ".sig_results.cpm", tag, ".csv", sep=""))
}

go_up_down <- function(genes, lengthData, GOmap, extra_tag = "UP", ...){
  cat(paste("\n\n### GO analysis with goseq - ", extra_tag, "\n", sep = ""))

  pwf <- nullp(genes, "hg38", "ensGene", bias.data = lengthData, plot.fit=FALSE)
  GO.wall <- goseq(pwf, gene2cat=GOmap)
  GO.wall$padj <- p.adjust(GO.wall$over_represented_pvalue, method="BH")
  GO.wall.sig <- subset(GO.wall, GO.wall$padj<.05)
  GO.wall.sig$ratio <- GO.wall.sig$numDEInCat / GO.wall.sig$numInCat
  
  if (nrow(GO.wall.sig) == 0) return()
  
  write.csv(GO.wall.sig, file=paste(directory, name, ".", cell, ".", contrast,
                                    ".sig_results.GO.", extra_tag, tag, ".csv", sep=""),
            row.names = FALSE)
  
  # Workaround in printing table in loop
  colSize <- ncol(GO.wall.sig)
  print(
    htmltools::tagList(
      datatable(GO.wall.sig, options = list(scrollX = TRUE), rownames = FALSE) %>%
        formatStyle(columns = c(1:colSize), fontSize = '80%')
    )
  )
  
  if (nrow(GO.wall.sig) >= 10){
    cat(paste(rep("\n", 11), collapse = ""))
  } else {
    cat(paste(rep("\n", nrow(GO.wall.sig)), collapse = ""))
  }
  
  cat("\n------------ Enriched", extra_tag, "-regulated GO ------------\r\n")
  # Print the details of the GO term
  enriched.GO <- GO.wall.sig$category
  print(enriched.GO)
  cat("\n\n")
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
  
  # Plot the GO terms in barplot, order the term by decresing numDEInCat
  GO.wall.sig$term <- reorder(GO.wall.sig$term, GO.wall.sig$numDEInCat)
  GO.wall.sig <- GO.wall.sig[order(GO.wall.sig$numDEInCat, decreasing=TRUE),]
  ggplot(subset(GO.wall.sig, !is.na(ontology)), 
         aes(term, numDEInCat, fill=ontology)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    coord_flip() + labs(title=paste(contrast, extra_tag, tag)) +
    scale_x_discrete(label=function(x) substr(x, 1, 40))

}

getPathwayName <- function(x) {
  x <- paste("hsa", x, sep = "")
  return(keggGet(x)[[1]]$NAME)
}

kegg_up_down <- function(genes, lengthData, KEGGmap, extra_tag = "UP", ...){
  cat(paste("\n\n\n### KEGG analysis with goseq - ", extra_tag, "\n", sep = ""))
  
  pwf <- nullp(genes, "hg38", "ensGene", bias.data = lengthData, plot.fit=FALSE)
  KEGG.wall <- goseq(pwf, gene2cat=KEGGmap)
  KEGG.wall$padj <- p.adjust(KEGG.wall$over_represented_pvalue, method="BH")
  KEGG.wall.sig <- subset(KEGG.wall, KEGG.wall$padj<.05)
  KEGG.wall.sig$ratio <- KEGG.wall.sig$numDEInCat / KEGG.wall.sig$numInCat
  
  if (nrow(KEGG.wall.sig) == 0) return()
  
  # Add pathway info to the table
  KEGG.wall.sig$term <- unlist(lapply(KEGG.wall.sig$category, getPathwayName))
  
  write.csv(KEGG.wall.sig, file=paste(directory, name, ".", cell, ".", contrast,
                                    ".sig_results.KEGG.", extra_tag, tag, ".csv", sep=""),
            row.names = FALSE)
  
  # Workaround in printing table in loop
  colSize <- ncol(KEGG.wall.sig)
  print(
    htmltools::tagList(
      datatable(KEGG.wall.sig, options = list(scrollX = TRUE), rownames = FALSE) %>%
        formatStyle(columns = c(1:colSize), fontSize = '80%')
    )
  )
  if (nrow(KEGG.wall.sig) >= 10){
    cat(paste(rep("\n", 11), collapse = ""))
  } else {
    cat(paste(rep("\n", nrow(KEGG.wall.sig)), collapse = ""))
  }
  
}

goseq_analysis <- function(...) {
  
  lengthData <- y$genes$length
  names(lengthData) <- rownames(y)
  genes <- as.vector(dt)
  names(genes) <- rownames(y)
  
  # UPs
  extra_tag <- "UP" 
  genes.up <- genes
  genes.up[genes.up == -1] <- 0 # Modify all the -1 to 0, look at the enrichment of the upregulated genes
  go_up_down(genes.up, lengthData, GOmap, extra_tag)
  kegg_up_down(genes.up, lengthData, KEGGmap, extra_tag)
  
  # DOWNs
  extra_tag <- "DOWN"
  genes.down <- genes
  genes.down[genes.down == 1] <- 0 
  genes.down[genes.down == -1] <- 1 # Modify all the -1 to 1, look at the enrichment of the downregulated genes
  go_up_down(genes.down, lengthData, GOmap, extra_tag)
  kegg_up_down(genes.down, lengthData, KEGGmap, extra_tag)
  

}



