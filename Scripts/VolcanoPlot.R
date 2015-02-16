library(ggplot2)
library(scales)
library(RColorBrewer)
require(grid)
library(DBI)
library(RMySQL)

#Functions to format plot and to scale axis
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
        theme(legend.key = element_rect(fill=color.background, colour=color.background)) +
        theme(legend.text = element_text(size=9,color=color.axis.title)) +
        theme(legend.title = element_blank()) +
        
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

#Get data from DB
con<-dbConnect(RMySQL::MySQL(), user='root', password='', host='localhost', dbname='SelaginellaGenomics')
Expression <- dbReadTable(conn = con,name = 'Expression')
DEgenes <- dbGetQuery(conn = con, statement = "SELECT * FROM Expression WHERE (rLogFC > 1 OR rLogFC < -1) AND DGEclust_padj < 0.01")

#need a couple of the formatting elements to draw lines
palette <- brewer.pal("Greys", n=9)
color.background = palette[2]
color.grid.major = palette[3]

#make the plot!
args <- commandArgs(trailingOnly = TRUE)
plotfile <- args[1]
png(plotfile, bg='transparent', res = 300, width=15, height=15, units='cm')
ggplot(Expression, aes(x=rLogFC, y=DGEclust_padj, colour=speciesID)) +
    geom_hline(yintercept=.01, colour="dark grey", size=1) +
    geom_segment(aes(x=-1, xend=1, y=.01, yend=.01), colour=color.background, size=1) +
    geom_segment(aes(x=-1, xend=1, y=.01, yend=.01), colour=color.grid.major, size=.1) +
    geom_segment(aes(x=-1, xend=-1, y=.01, yend=.00001), colour="dark grey") +
    geom_segment(aes(x=1, xend=1, y=.01, yend=.00001), colour="dark grey") +
    geom_point(alpha=0.2, size=1) +
    geom_point(data = DEgenes, size=2, aes(colour=speciesID)) +
    scale_y_continuous(trans = reverselog_trans(10), breaks=c(.1,.01,.001, 0.0001, 0.00001), labels=comma) +
    scale_x_continuous(breaks=c(-3, -2, -1, 0, 1, 2, 3)) +
    labs(x="fold change (log2)", y="p-values") +
    fte_theme() +
    theme(legend.position="right") +
    guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()
