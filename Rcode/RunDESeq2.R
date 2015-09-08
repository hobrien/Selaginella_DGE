library(DESeq2)


args <- commandArgs(trailingOnly = TRUE)
CountTable =  read.table( file=args[1] , header=TRUE, row.names=1 )
sample1 <- args[2]
sample2 <- args[3]
textfile <- args[4] 
species <- substr(sample1, 1, nchar(sample1) - 1)
name = paste(sample1, substr(sample2, nchar(sample2), nchar(sample2)), sep='')
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
write.table(results, textfile)
