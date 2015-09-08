Selaginella_DGE
===============

Differential Gene Expression Analyses on Selaginella

The Script SelaginellaPipeline.sh has all (or probably most) of the commands to do the following:

- Run Transdecoder to extract coding sequences from Trinity contigs
- Run an all-by-all blastp analysis of amino acid translations
- Map reads to Trinity contigs keeping all alignments
- Run Corset to cluster contigs
- Create MySQL database with contigs, blast results and Corset results
- Extract sequence from each Corset cluster with highest blast hit and write to file
- Run OrthoMCL on highest scoring sequences and add results to database
- Use count data from Corset to run differential expression analyses using DESeq and DGEclust
- Make trees for each gene family and map on expression data

Dependancies
------------
- Blast+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- Corset (https://code.google.com/p/corset-project/)
- DESeq (http://bioconductor.org/packages/release/bioc/html/DESeq.html)
- DGEClust (http://dvav.me/dgeclust/)
- ETE2 (http://etetoolkit.org/)
- MAFFT (http://mafft.cbrc.jp/alignment/software/)
- MySQL (http://www.mysql.com/)
- OrthoMCL (http://www.orthomcl.org/orthomcl/)
- PhyML (http://atgc.lirmm.fr/phyml/)
- R (http://www.r-project.org/)
- Transdecoder (http://transdecoder.sourceforge.net/)
- Trinity (http://trinityrnaseq.sourceforge.net/)

Files
-----
Adaptor sequence and low quality base pairs were trimmed using the fastq-mcf utility from ea-utils {eautilsCommand:6tDggvj9} using a quality threshold of 10.


