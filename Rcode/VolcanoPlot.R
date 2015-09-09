library(ggplot2)
library(scales)
library(RColorBrewer)
require(grid)
library(DBI)
library(RMySQL)

source("~/Documents/BTsync2/R/FormatGGplot.R")
#Get data from DB
con<-dbConnect(RMySQL::MySQL(), user='root', password='', host='localhost', dbname='SelaginellaGenomics')
Expression <- dbReadTable(conn = con,name = 'Expression')
#Expression <- dbGetQuery(conn = con, statement = "SELECT * FROM Expression WHERE sample1 = 'KRAUS1' AND sample2 = 'KRAUS2'")
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
    geom_hline(yintercept=.01, colour="dark grey") +
    geom_segment(aes(x=-1, xend=1, y=.01, yend=.01), colour=color.background, size=1) +
    geom_segment(aes(x=-1, xend=1, y=.01, yend=.01), colour=color.grid.major, size=.1) +
    geom_segment(aes(x=-1, xend=-1, y=.01, yend=.00001), colour="dark grey") +
    geom_segment(aes(x=1, xend=1, y=.01, yend=.00001), colour="dark grey") +
    geom_vline(xintercept=0, colour=color.grid.major, size=.1) +
    geom_point(alpha=0.2, size=1) +
    geom_point(data = DEgenes, size=1, aes(colour=speciesID)) +
    scale_y_continuous(trans = reverselog_trans(10), breaks=c(.1,.01,.001, 0.0001, 0.00001), labels=comma) +
    scale_x_continuous(breaks=c(-3, -2, -1, 0, 1, 2, 3)) +
    labs(x="fold change (log2)", y="p-values") +
    fte_theme() +
    theme(legend.position="right") +
    guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()
