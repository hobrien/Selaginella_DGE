library(ggplot2)
library(plyr)
library(DESeq2)
library(scales)
library(RColorBrewer)
require(grid)

fte_theme <- function() {
    
    # Generate the colors for the chart procedurally with RColorBrewer
    palette <- brewer.pal("Greys", n=9)
    color.background = palette[2]
    color.grid.major = palette[3]
    color.axis.text = palette[6]
    color.axis.title = palette[7]
    color.title = palette[9]
    
    # Begin construction of chart
    theme_bw(base_size=9) +
        
        # Set the entire chart region to a light gray color
        theme(panel.background=element_rect(fill=color.background, color=color.background)) +
        theme(plot.background=element_rect(fill=color.background, color=color.background)) +
        theme(panel.border=element_rect(color=color.background)) +
        
        # Format the grid
        theme(panel.grid.major=element_line(color=color.grid.major,size=.25)) +
        theme(panel.grid.minor=element_blank()) +
        theme(axis.ticks=element_blank()) +
        
        # Format the legend, but hide by default
        theme(legend.position="none") +
        theme(legend.background = element_rect(fill=color.background)) +
        theme(legend.text = element_text(size=7,color=color.axis.title)) +
        
        # Set title and axis labels, and format these and tick marks
        theme(plot.title=element_text(color=color.title, size=10, vjust=1.25)) +
        theme(axis.text.x=element_text(size=9,color=color.axis.text)) +
        theme(axis.text.y=element_text(size=9,color=color.axis.text)) +
        theme(axis.title.x=element_text(size=10,color=color.axis.title, vjust=0)) +
        theme(axis.title.y=element_text(size=10,color=color.axis.title, vjust=1.25)) +
        
        # Plot margins
        theme(plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm"))
}
reverselog_trans <- function(base = exp(1)) {
     trans <- function(x) -log(x, base)
     inv <- function(x) base^(-x)
     trans_new(paste0("reverselog-", format(base)), trans, inv, 
               log_breaks(base = base), 
               domain = c(1e-100, Inf))
}


args <- commandArgs(trailingOnly = TRUE)
DGEclust_file <- args[1] #"Bioinformatics/Selaginella_DGE/DGEclust/KRAUS_12.txt"
#CountTable =  read.table( file="Bioinformatics/Selaginella_DGE/DGEclust/KRAUScounts.txt" , header=TRUE, row.names=1 )
CountTable =  read.table( file=args[2] , header=TRUE, row.names=1 )
#DGEclust_cutoff <- as.numeric(0.01)
DGEclust_cutoff <- as.numeric(args[3])
plotfile<-paste(strsplit(DGEclust_file, '\\.')[[1]][1], 'volcano.pdf', sep='_') 
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
rld <- rlogTransformation( DEseq_data )
results <- data.frame(
    assay(rld), 
    avgLogExpr = ( assay(rld)[,2] + assay(rld)[,1] ) / 2,
    rLogFC = assay(rld)[,2] - assay(rld)[,1] )
DGEclust<-read.table(file=DGEclust_file, header=T)
results$"DGEclust_padj"<-DGEclust[,"padj"][match(row.names(results), row.names(DGEclust))]

write.table(results[results$DGEclust_padj < DGEclust_cutoff & abs(results$rLogFC)>1, ], textfile)
title <- sprintf("Volcano plot for %s vs %s", sample1, sample2)
pdf(plotfile)
sig_res<-results[results$DGEclust_padj < 0.01 & abs(results$rLogFC)>1,]
if (nrow(sig_res) > 0){
	ggplot(results, aes(x=rLogFC, y=DGEclust_padj)) +
		geom_point(alpha=0.1, size=2, colour="dark blue") +
		geom_point(size=2, data=sig_res, aes(colour='red')) +
		fte_theme() +
		scale_y_continuous(trans = reverselog_trans(10), breaks=c(.1,.01,.001, 0.0001), labels=comma) +
		labs(title=title, x="fold change (log2)", y="p-values")
	} else {
	ggplot(results, aes(x=rLogFC, y=DGEclust_padj)) +
		geom_point(alpha=0.1, size=2, colour="dark blue") +
		fte_theme() +
		scale_y_continuous(trans = reverselog_trans(10), breaks=c(.1,.01,.001, 0.0001), labels=comma) +
		labs(title=title, x="fold change (log2)", y="p-values")
}
dev.off()


