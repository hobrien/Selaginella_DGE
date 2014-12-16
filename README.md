Selaginella_DGE
===============

Differential Gene Expression Analyses on Selaginella

The Script SelaginellaPipeline.sh has all (or probably most) of the commands to do the following:

-Run Transdecoder to extract coding sequences from Trinity contigs
-Run an all-by-all blastp analysis of amino acid translations
-Map reads to Trinity contigs keeping all alignments
-Run Corset to cluster contigs
-Create MySQL database with contigs, blast results and Corset results
-Extract sequence from each Corset cluster with highest blast hit and write to file
-Run OrthoMCL on highest scoring sequences and add results to database
-Use count data from Corset to run differential expression analyses using DESeq and DGEclust
-Make trees for each gene family and map on expression data

Dependancies
------------
-Blast+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
-Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
-Corset (https://code.google.com/p/corset-project/)
-DESeq (http://bioconductor.org/packages/release/bioc/html/DESeq.html)
-ETE2 (http://etetoolkit.org/)
-MAFFT (http://mafft.cbrc.jp/alignment/software/)
-MySQL (http://www.mysql.com/)
-OrthoMCL (http://www.orthomcl.org/orthomcl/)
-PhyML (http://atgc.lirmm.fr/phyml/)
-R (http://www.r-project.org/)
-Transdecoder (http://transdecoder.sourceforge.net/)
-Trinity (http://trinityrnaseq.sourceforge.net/)

Files
-----

All_by_all:
clust -t 10000 -dt 10 -o All_by_All Selag_counts.txt 

By_leaf:
clust -t 10000 -dt 10 -o By_leaf -g [[0,4,7,11],[1,5,8,12],[2,6,9,13],[3,10,14]] Selag_counts.txt

By_species:
clust -t 10000 -dt 10 -o By_species -g [[0,1,2,3],[4,5,6],[7,8,9,10],[11,12,13,14]] Selag_counts.txt 

GetPvals.py: script to calculate p vlues for all comparisons

Selag_counts.txt: count data for expression analyses


