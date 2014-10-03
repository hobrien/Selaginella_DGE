library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
infilename <- args[1]
outfilename <- args[2]
infile <- read.table(infilename, sep='\t', header=F)
trans_theme <- theme(
  panel.grid.major = element_line(color='black'),
  panel.grid.minor = element_line(color=NA),
  panel.background = element_rect(fill=NA, colour='black'),
  plot.background = element_rect(fill=NA),
  axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
  axis.text.y = element_text(colour="black",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
  axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
  axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),
  strip.text = element_text(colour=c("black", "green", "red", "orange", "blue"), size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
  strip.background = element_rect(fill=NA)
)

png(file=outfilename, bg = "transparent", width=1000, height=1000)
ggplot(infile, aes(V3, fill =V1)) + 
  geom_histogram(binwidth=.1) +
  facet_wrap(~ V1) +
  scale_x_log10(name="length") +
  scale_fill_manual(values=c('green', 'red', 'orange', 'blue')) +
  theme(legend.position="none") +
  trans_theme
dev.off()
