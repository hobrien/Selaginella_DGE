library(ggplot2)
library(plyr)
library(DESeq2)
library(VennDiagram)


MergeResults<-function(DEseq, DGEclust, name){
  #rename columns and add info about species and comparison
  #row.names(DEseq)<-DEseq$id
  DEseq["samples"]<-name
  DEseq<-rename(DEseq, c("pvalue" = "DEseq_pval", "padj"="DEseq_padj"))
  
  #rename DGEclust columns, if necessary
  #DGEclust<-rename(DGEclust, c("Posteriors" = "pval", "FDR" = "padj"))
  DEseq$"DGEclust_pval"<-DGEclust[,"pval"][match(row.names(DEseq), row.names(DGEclust))]
  DEseq$"DGEclust_padj"<-DGEclust[,"padj"][match(row.names(DEseq), row.names(DGEclust))]
  
  
  #remove missing values that will cause an error in the following step
  DEseq<-na.omit(DEseq)
  
  #this will replace extreme values with zero, which makes the plot nicer
  if (nrow(DEseq[DEseq$DEseq_padj < 1e-5, ]["DEseq_padj"]) > 0 ) {
    DEseq[DEseq$DEseq_padj < 1e-5, ]["DEseq_padj"] <- 0
  }
  return(as.data.frame(DEseq))
}
scatterplot<-function(pvalues, DEseq_cutoff, DGEclust_cutoff, name) {
  #make scatterplot of adjusted p-values and Venn Diagram of overlap
  sp = ggplot(pvalues, aes(x=DEseq_padj, y=DGEclust_padj)) +
    geom_point(size=1.2) +
    scale_x_log10(limits=c(1e-5,1), breaks=c(1e-5, 1e-4, 1e-3, 0.01, 0.1, 1)) +
    scale_y_log10(limits=c(1e-5,1), breaks=c(1e-5, 1e-4, 1e-3, 0.01, 0.1, 1)) +
    geom_hline(yintercept=DGEclust_cutoff, size=.2) +
    geom_vline(xintercept=DEseq_cutoff, size=.2) +
    annotate("text", x=1e-5, y=1, hjust=0, label=sprintf("r = %.3f", cor(results$DEseq_padj,results$DGEclust_padj))) +
    annotate("text", hjust=0, x=1e-5, y=.7, label=sprintf("p = %.5f", cor.test(pvalues$DEseq_padj,pvalues$DGEclust_padj)$p.value)) +
    theme_bw() +
    ggtitle(name)
  #geom_point(data=results[results$id == 'cluster_19641',], colour='red')
  return(sp)
}
vennplot <-function(pvalues, DEseq_cutoff, DGEclust_cutoff) {
  if ( nrow(pvalues[pvalues$DGEclust_padj < DGEclust_cutoff, ]) > 0 && nrow(pvalues[pvalues$DEseq_padj < DEseq_cutoff, ]) > 0 ) {
    if ( nrow(pvalues[pvalues$DGEclust_padj < DGEclust_cutoff, ]) < nrow(pvalues[pvalues$DEseq_padj < DEseq_cutoff, ]) ) {
      draw.pairwise.venn(area1=nrow(pvalues[results$DGEclust_padj < DGEclust_cutoff, ]), 
                         area2=nrow(pvalues[pvalues$DEseq_padj < DEseq_cutoff, ]), 
                         cross.area=nrow(pvalues[pvalues$DEseq_padj < DEseq_cutoff & pvalues$DGEclust_padj < DGEclust_cutoff, ]), 
                         fill=c('blue', 'red'), 
                         cat.col = c('blue', 'red'),
                         category=c("DEseq", "DGEclust"), 
                         lty = 'blank',
                         cex = 2,
                         cat.cex = 1.75,
                         margin=0.2,
                         fontfamily='sans',
                         cat.fontfamily='sans')
    }
    else {
      draw.pairwise.venn(area1=nrow(pvalues[pvalues$DGEclust_padj < DGEclust_cutoff, ]),
                         area2=nrow(pvalues[pvalues$DEseq_padj < DEseq_cutoff, ]), 
                         cross.area=nrow(pvalues[pvalues$DEseq_padj < DEseq_cutoff & pvalues$DGEclust_padj < DGEclust_cutoff, ]), 
                         fill=c('red', 'blue'), 
                         cat.col = c('red', 'blue'),
                         category=c("DGEclust", "DEseq"), 
                         lty = 'blank',
                         cex = 2,
                         cat.cex = 1.75,
                         margin=0.2,
                         fontfamily='sans',
                         cat.fontfamily='sans')
    }
  }
  else if (nrow(pvalues[pvalues$DGEclust_padj < DGEclust_cutoff, ]) > 0) {
    draw.single.venn(area=nrow(pvalues[pvalues$DGEclust_padj < DGEclust_cutoff, ]), 
                     fill=c('red'), 
                     cat.col = c('red'),
                     category=c("DGEclust"), 
                     lty = 'blank',
                     cex = 2,
                     cat.cex = 1.75,
                     margin=0.2,
                     fontfamily='sans',
                     cat.fontfamily='sans')
  }
  else if (nrow(pvalues[pvalues$DEseq_padj < DEseq_cutoff, ]) > 0) {
    draw.single.venn(area=nrow(pvalues[pvalues$DEseq_padj < DEseq_cutoff, ]), 
                     fill=c('blue'), 
                     cat.col = c('blue'),
                     category=c("DEseq"), 
                     lty = 'blank',
                     cex = 2,
                     cat.cex = 1.75,
                     margin=0.2,
                     fontfamily='sans',
                     cat.fontfamily='sans')
  }
}

args <- commandArgs(trailingOnly = TRUE)
DGEclust_file <- args[1] #"Bioinformatics/Selaginella_DGE/DGEclust/KRAUS_12.txt"
#CountTable =  read.table( file="Bioinformatics/Selaginella_DGE/DGEclust/KRAUScounts.txt" , header=TRUE, row.names=1 )
CountTable =  read.table( file=args[2] , header=TRUE, row.names=1 )
#DGEclust_cutoff <- as.numeric(0.01)
DGEclust_cutoff <- as.numeric(args[3])
DEseq_cutoff<- DGEclust_cutoff
plotfile<-paste(strsplit(DGEclust_file, '\\.')[[1]][1], DGEclust_cutoff, 'pvals2.pdf', sep='_') 
textfile<-paste(strsplit(DGEclust_file, '\\.')[[1]][1], DGEclust_cutoff, 'overlap2.txt', sep='_') 
species_regex<-gregexpr('(KRAUS)|(MOEL)|(UNC)|(WILD)', DGEclust_file)
regex_start <- species_regex[[1]][1]
regex_end <- regex_start + attr(species_regex[[1]], 'match.length')[1] -1
species <- substr(DGEclust_file, regex_start, regex_end)
sample1 <- paste(species, substr(DGEclust_file, regex_end + 2, regex_end + 2), sep='')
sample2 <- paste(species, substr(DGEclust_file, regex_end + 3, regex_end + 3), sep='')
name = substr(DGEclust_file, regex_start, regex_end + 3)
CountTable<-CountTable[,c(sample1, sample2)]
Setup<- data.frame(row.names=colnames(CountTable), condition = factor( c( sample1, sample2 )))
DEseq_data<-DESeqDataSetFromMatrix(countData = CountTable,
                                  colData = Setup,
                                  design = ~ condition)
DEseq_data <- DESeq(DEseq_data)
DEseq_results <- results(DEseq_data)
DGEclust<-read.table(file=DGEclust_file, header=T)
results<-MergeResults(DEseq_results, DGEclust, name)
#write list of overlaps to textfile
write.table(results[results$DEseq_padj < DEseq_cutoff & results$DGEclust_padj < DGEclust_cutoff, ], textfile)

pdf(plotfile)
print(scatterplot(results, DEseq_cutoff, DGEclust_cutoff, name))
print(vennplot(results, DEseq_cutoff, DGEclust_cutoff))
dev.off()


