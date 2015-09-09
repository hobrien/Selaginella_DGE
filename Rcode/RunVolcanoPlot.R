source("/Users/HeathOBrien/Bioinformatics/Selaginella_DGE/Rcode/VolcanoPlot.R")
args <- commandArgs(trailingOnly = TRUE)
plotfile <- args[1]
png(plotfile, bg='transparent', res = 300, width=15, height=15, units='cm')
print(Volcanoplot(ALL))
dev.off()