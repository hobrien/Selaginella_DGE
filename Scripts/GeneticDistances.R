library("ape")
library("reshape2")
library("ggplot2")

GetDist <- function(distance, seq1, seq2){
  selection<-paste("'", seq1, "'", ", ", "'", seq2, "'", sep="")
  eval(parse(text = paste("as.matrix(distance)[",
                          selection, "]")))
}

setwd("~/Bioinformatics/Selaginella/Alignments")
Species = c("KRAUS", "MOEL", "UNC", "WILD")
results<-data.frame(Comparison=character(), Distance=numeric(), Locus=character())
for (file in list.files('.')) {
  print(file) 
  aln<-read.dna(file=file, format="fasta")
  dist<-dist.dna(aln, model="raw", pairwise.deletion=T)
  attr(dist, "Labels")<-as.list(colsplit(string=attr(dist, "Labels"), pattern="\\Q|\\E", names=c(1,2))[1])[[1]]
  for (x in 1:3){
    for (y in (x+1):4){
      if (Species[x] %in% labels(dist) && Species[y] %in% labels(dist)) {
        results<-rbind(results,data.frame(Comparison=paste(Species[x], Species[y], sep="-"), Distance=GetDist(dist, Species[x], Species[y]), Locus=file))
      }
    }  
  }
}
pdf("~/Google Drive/Selaginella/Results/Pairwise_distances.pdf")
ggplot(results, aes(Distance))+geom_histogram()+facet_wrap(~Comparison)+theme_bw()
dev.off()
