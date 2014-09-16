library(ggplot2)
library(plyr)
library(DESeq)
library(VennDiagram)


MergeResults<-function(DEseq, DGEclust, name, textfile){
  #rename columns and add info about species and comparison
  DEseq["samples"]<-name
  DEseq<-rename(DEseq, c("pval" = "DESeq_pval", "padj"="DESeq_padj"))

  #sort rows
  DGEClust<-DGEClust[order(as.numeric(rownames(DGEClust))),]
  
  #rename DGEClust columns, if necessary
  DGEClust<-rename(DGEClust, c("Posteriors" = "pval", "FDR" = "padj"))
  #add DGESClust data to results
  DEseq["DGEClust_pval"]<-DGEClust["pval"]
  DEseq["DGEClust_padj"]<-DGEClust["padj"]

  #write list of overlaps to textfile
  write.table(DEseq[DEseq$DESeq_padj < 0.1 & DEseq$DGEClust_padj < 0.1, ], textfile)

  #remove missing values that will cause an error in the following step
  DEseq<-na.omit(DEseq)

  #this will replace extreme values with zero, which makes the plot nicer
  if (nrow(DEseq[DEseq$DESeq_padj < 1e-5, ]["DESeq_padj"]) > 0 ) {
    DEseq[DEseq$DESeq_padj < 1e-5, ]["DESeq_padj"] <- 0
  }
  return(DEseq)
}
scatterplot<-function(pvalues, DESeq_cutoff, DGEClust_cutoff, name) {
  #make scatterplot of adjusted p-values and Venn Diagram of overlap
  sp = ggplot(pvalues, aes(x=DESeq_padj, y=DGEClust_padj)) +
    geom_point(size=1.2) +
    scale_x_log10(limits=c(1e-5,1), breaks=c(1e-5, 1e-4, 1e-3, 0.01, 0.1, 1)) +
    scale_y_log10(limits=c(1e-5,1), breaks=c(1e-5, 1e-4, 1e-3, 0.01, 0.1, 1)) +
    geom_hline(yintercept=DGEClust_cutoff, size=.2) +
    geom_vline(xintercept=DESeq_cutoff, size=.2) +
    annotate("text", x=1e-5, y=1, hjust=0, label=sprintf("r = %.3f", cor(results$DESeq_padj,results$DGEClust_padj))) +
    annotate("text", hjust=0, x=1e-5, y=.7, label=sprintf("p = %.5f", cor.test(pvalues$DESeq_padj,pvalues$DGEClust_padj)$p.value)) +
    theme_bw() +
    ggtitle(name)
    #geom_point(data=results[results$id == 'cluster_19641',], colour='red')
    return(sp)
}
vennplot <-function(pvalues, DESeq_cutoff, DGEClust_cutoff) {
   if ( nrow(pvalues[pvalues$DGEClust_padj < DGEClust_cutoff, ]) > 0 && nrow(pvalues[pvalues$DESeq_padj < DESeq_cutoff, ]) > 0 ) {
      if ( nrow(pvalues[pvalues$DGEClust_padj < DGEClust_cutoff, ]) < nrow(pvalues[pvalues$DESeq_padj < DESeq_cutoff, ]) ) {
         draw.pairwise.venn(area1=nrow(pvalues[results$DGEClust_padj < DGEClust_cutoff, ]), 
            area2=nrow(pvalues[pvalues$DESeq_padj < DESeq_cutoff, ]), 
            cross.area=nrow(pvalues[pvalues$DESeq_padj < DESeq_cutoff & pvalues$DGEClust_padj < DGEClust_cutoff, ]), 
            fill=c('blue', 'red'), 
            cat.col = c('blue', 'red'),
            category=c("DESeq", "DGEClust"), 
            lty = 'blank',
            cex = 2,
            cat.cex = 1.75,
            margin=0.2,
            fontfamily='sans',
            cat.fontfamily='sans')
      }
      else {
         draw.pairwise.venn(area1=nrow(pvalues[pvalues$DESeq_padj < DESeq_cutoff, ]), 
            area2=nrow(pvalues[pvalues$DGEClust_padj < DGEClust_cutoff, ]), 
            cross.area=nrow(pvalues[pvalues$DESeq_padj < DESeq_cutoff & pvalues$DGEClust_padj < DGEClust_cutoff, ]), 
            fill=c('red', 'blue'), 
            cat.col = c('red', 'blue'),
            category=c("DGEClust", "DESeq"), 
            lty = 'blank',
            cex = 2,
            cat.cex = 1.75,
            margin=0.2,
            fontfamily='sans',
            cat.fontfamily='sans')
      }
   }
   else if (nrow(pvalues[pvalues$DGEClust_padj < DGEClust_cutoff, ]) > 0) {
      draw.single.venn(area=nrow(pvalues[pvalues$DGEClust_padj < DGEClust_cutoff, ]), 
                       fill=c('red'), 
                       cat.col = c('red'),
                       category=c("DGEClust"), 
                       lty = 'blank',
                       cex = 2,
                       cat.cex = 1.75,
                       margin=0.2,
                       fontfamily='sans',
                       cat.fontfamily='sans')
   }
   else if (nrow(pvalues[pvalues$DESeq_padj < DESeq_cutoff, ]) > 0) {
      draw.single.venn(area=nrow(pvalues[pvalues$DESeq_padj < DESeq_cutoff, ]), 
                       fill=c('blue'), 
                       cat.col = c('blue'),
                       category=c("DESeq"), 
                       lty = 'blank',
                       cex = 2,
                       cat.cex = 1.75,
                       margin=0.2,
                       fontfamily='sans',
                       cat.fontfamily='sans')
   }
}
 
#read in raw count data
selaginellaCountTable =  read.table( file="/Users/HeathOBrien/Bioinformatics/Selaginella/Counts/Selag_counts.txt" , header=TRUE, row.names=1 )
species.list <- c('KRAUS', 'MOEL', 'UNC', 'WILD')

#This is going to be a major overhall to make it work with the specified file.
#Can be looped using a shell script


args <- commandArgs(trailingOnly = TRUE)
DGEClust_file <-args[1]
DGEClust_cutoff <- as.numeric(args[2])
DESeq_cutoff<- DGEClust_cutoff
plotfile<-paste(strsplit(DGEClust_file, '\\.')[[1]][1], DGEClust_cutoff, 'pdf', sep='.') 
textfile<-paste(strsplit(DGEClust_file, '\\.')[[1]][1], DGEClust_cutoff, 'overlap.txt', sep='_') 
species_regex<-gregexpr('(KRAUS)|(MOEL)|(UNC)|(WILD)', DGEClust_file)
if ( length(species_regex[[1]]) == 2 ) {  #Analyse Species comparisons
   sample1 <- substr(DGEClust_file, species_regex[[1]][1], species_regex[[1]][1] + attr(species_regex[[1]], 'match.length')[1] -1)
   sample2 <- substr(DGEClust_file, species_regex[[1]][2], species_regex[[1]][2] + attr(species_regex[[1]], 'match.length')[2] -1)
   name = paste(sample1, sample2, sep="_")
   condition = factor( c( 'KRAUS', 'KRAUS', 'KRAUS', 'KRAUS', 'MOEL', 'MOEL', 'MOEL', 'UNC', 'UNC', 'UNC', 'UNC', 'WILD', 'WILD', 'WILD', 'WILD' ) )
   DESeq<-newCountDataSet( selaginellaCountTable, condition)
   DESeq <- estimateSizeFactors(DESeq)
   DESeq <- estimateDispersions(DESeq)
} else if ( species_regex[[1]][1] > -1 ) { #Analyse comparison between two leaves within a species
   regex_start <- species_regex[[1]][1]
   regex_end <- regex_start + attr(species_regex[[1]], 'match.length')[1] -1
   species <- substr(DGEClust_file, regex_start, regex_end)
   sample1 <- paste("leaf", substr(DGEClust_file, regex_end + 1, regex_end + 1), sep='')
   sample2 <- paste("leaf", substr(DGEClust_file, regex_end + 2, regex_end + 2), sep='')
   name = substr(DGEClust_file, regex_start, regex_end + 2)
   if (species == 'MOEL') {
     condition = factor( c( "leaf1", "leaf2", "leaf3" ) )
     CountTable = selaginellaCountTable [c(paste(species, 1, sep=''), paste(species, 2, sep=''), paste(species, 3, sep=''))]
  }
  else {
     condition = factor( c( "leaf1", "leaf2", "leaf3", "leaf4" ) )
     CountTable = selaginellaCountTable[c(paste(species, 1, sep=''), paste(species, 2, sep=''), paste(species, 3, sep=''), paste(species, 4, sep=''))]
  }
  DESeq<-newCountDataSet( CountTable, condition)
  DESeq <- estimateSizeFactors(DESeq)
  DESeq <- estimateDispersions(DESeq, method="blind", sharingMode="fit-only")
} else {
   leaf_regex<-gregexpr('(leaf)|(all)', DGEClust_file, ignore.case=T)
   if ( leaf_regex[[1]][1] > -1 ) { #Analyse comparisons between leaves pooled by species
      regex_start <- leaf_regex[[1]][1]
      regex_end <- regex_start + attr(leaf_regex[[1]], 'match.length')[1] -1
      species <- 'ALL'
      sample1 <- paste("leaf", substr(DGEClust_file, regex_end + 1, regex_end + 1), sep='')
      sample2 <- paste("leaf", substr(DGEClust_file, regex_end + 2, regex_end + 2), sep='')
      name = substr(DGEClust_file, regex_start, regex_end + 2)
      condition = factor( c( "leaf1", "leaf2", "leaf3", "leaf4", "leaf1", "leaf2", "leaf3", "leaf1", "leaf2", "leaf3", "leaf4", "leaf1", "leaf2", "leaf3", "leaf4") )
      DESeq<-newCountDataSet(selaginellaCountTable, condition)
      DESeq <- estimateSizeFactors(DESeq)
      DESeq <- estimateDispersions(DESeq)
  }
}
DESeq.results <- nbinomTest(DESeq, sample1, sample2)
DGEClust<-read.table(file=DGEClust_file, header=T)
DGEClust<-rename(DGEClust, c("Posteriors" = "pval", "FDR" = "padj"))
results<-MergeResults(DESeq.results, DGEClust, name, textfile)
pdf(plotfile)
print(scatterplot(results, DESeq_cutoff, DGEClust_cutoff, name))
print(vennplot(results, DESeq_cutoff, DGEClust_cutoff))
dev.off()
       
    
