MAXCLUST=30466
BASEDIR=/Users/HeathOBrien/Bioinformatics

#initialise database with species info
mysql -u root SelaginellaGenomics < /Users/HeathOBrien/Google\ Drive/Selaginella/DB/selaginellaSQL.txt

#Add BLUELEAF sequences to database
for file in `find /Users/HeathOBrien/Bioinformatics/Selaginella/Assemblies | grep "_Tr.fa$"`
do
	AddData.py -f sequences -i $file
done

#Add reference transcriptome sequences to database
for file in `find /Users/HeathOBrien/Bioinformatics/Selaginella/RefSeq -name *_cds.fa`
do
	AddData.py -f sequences -i $file
done

AddData.py -f sequences -i /Users/HeathOBrien/Bioinformatics/Selaginella/RefSeq/Selaginella_moellendorffii.v1.0.17.cdna.all.fa

#Add info about coding regions of transcripts to database
for file in `find /Users/HeathOBrien/Bioinformatics/Selaginella/Transdecoder | grep _Tr.fa.transdecoder.bed`
do
	AddData.py -f coding -i $file
done

AddData.py -f ref_coding #reference sequences are all coding, so we don't need info from BED file

#Add orthoMCL info to DB
AddData.py -f orthologs -i /Users/HeathOBrien/Bioinformatics/Selaginella/OrthoMCL/Selaginella_groups.txt

######################################### RUN PHYLOTREE PRUNER #############################################
mkdir /Users/HeathOBrien/Bioinformatics/Selaginella/Clusters
cd /Users/HeathOBrien/Bioinformatics/Selaginella/Clusters

#Write sequences from each ortholog group to files
for i in $(seq 0 $MAXCLUST)
do
  GetData.py -f seq_clusters -c $i > cluster.$i.fa
done

#run translatorX on all fasta files
find . -name "cluster.*.fa"| xargs -P 6 -n 1 translatorx_vLocal.pl -p F -i

#convert to Phylip (this part isn't parallelised, but it shouldn't take too long)
for i in `ls | grep nt_ali.fasta$| cut -d. -f2`
do
  trimal -gappyout -in cluster.$i.nt_ali.fasta | ConvertAln.py -x fasta -f phylip -i STDIN -o cluster.$i.phy
done

#run phyml on all phylip files
find . -name "cluster.*.phy"| xargs -P 6 -n 1 phyml --quiet --no_memory_check -o n -b 0 -i

#run phylotree pruner on results and add results to DB
for i in `ls | grep .phy$| cut -d. -f2`
do
  mv cluster.$i.phy_phyml_tree.txt cluster.$i.nwk
  java PhyloTreePruner cluster.$i.nwk 4 cluster.$i.fa 0.5 u
  AddData.py -f nr -i cluster.$i.fa_pruned.fa
done

#clean up unwanted files 
find . -name "*_ali.fasta" |xargs rm
find . -name "*.html" |xargs rm
find . -name "*aaseqs*" |xargs rm
find . -name "*.log" |xargs rm
find . -name "*_pruned.fa" | xargs rm
find . -name "*_phyml_stats.txt" | xargs rm

####################################### MAKE TREES FOR EACH NR SET ###########################################
mkdir /Users/HeathOBrien/Bioinformatics/Selaginella/NRclusters
cd /Users/HeathOBrien/Bioinformatics/Selaginella/NRclusters

#Write sequences from each ortholog group to files
for i in $(seq 0 $MAXCLUST)
do
  GetData.py -f nr_clusters -c $i -o cluster.$i.fa
done

#run translatorX on all fasta files
find . -name "cluster.*.fa"| xargs -P 6 -n 1 translatorx_vLocal.pl -p F -i

#convert to Phylip (this part isn't parallelised, but it shouldn't take too long)
for i in `ls | grep nt_ali.fasta$| cut -d. -f2`
do
  trimal -gappyout -in cluster.$i.nt_ali.fasta | ConvertAln.py -x fasta -f phylip -i STDIN -o cluster.$i.phy
done

#run phyml on all phylip files
find . -name "cluster.*.phy"| xargs -P 6 -n 1 phyml --quiet --no_memory_check -o n -b 0 -i

#rename tree files
for i in `ls | grep .phy$| cut -d. -f2`
do
  mv cluster.$i.phy_phyml_tree.txt cluster.$i.nwk
done

#clean up unwanted files 
find . -name "*_ali.fasta" |xargs rm
find . -name "*.html" |xargs rm
find . -name "*aaseqs*" |xargs rm
find . -name "*.log" |xargs rm
find . -name "*_pruned.fa" | xargs rm
find . -name "*_phyml_stats.txt" | xargs rm

####################################### MAKE GRAPH OF GENETIC DISTANCES ###########################################
Counts=~/Bioinformatics/Selaginella/Counts
mkdir /Users/HeathOBrien/Bioinformatics/Selaginella/Alignments
cd /Users/HeathOBrien/Bioinformatics/Selaginella/Alignments

for i in `ls ../NRclusters | grep .phy | cut -d. -f2`
do 
  if [[ -s cluster.$i.phy ]]
  then
    ConvertAln.py -i ../NRclusters/cluster.$i.phy -o cluster.$i.fa -f fasta; done;
  else
    rm cluster.$i.phy
    rm cluster.$i.nwk
  fi
done

Rscript ~/Documents/R/GeneticDistances.R

#################################### MULTIMAP READS AND RUN CORSET ####################################
for species in 'KRAUS' 'MOEL' 'UNC' 'WILD'
do
  multimap.sh
  RunCorset.sh
done

#################################### UPLOAD BLAST DATA ####################################
for species in 'KRAUS' 'MOEL' 'UNC' 'WILD'
do
  if test $species == 'KRAUS'
  then
    grep ^KRUS $BASEDIR/OrthoMCL/goodProteins.bl > $species
  else
    grep ^$species $BASEDIR/OrthoMCL/goodProteins.bl > $species
  fi
  $BASEDIR/Scripts/AddData.py -f blast -i species
done

################################ MAP READS TO NR SEQS AND GET HIT COUNTS ##################################
cd /Users/HeathOBrien/Bioinformatics/Selaginella/Non_redundant
for species in 'KRAUS' 'MOEL' 'UNC' 'WILD'
do
  GetData.py -f nr_seqs -s $species > ${species}_nr.fa
  GetData.py -f nr_seq -s $species > ${species}
done
CompareSets.pl UNC MOEL KRAUS WILD
mv VennPlot.pdf ~/Google\ Drive/Selaginella/Figures/OrthoVenn.pdf
  
bash ~/Bioinformatics/Selaginella/Shell_scripts/unimap.sh KRAUS MOEL UNC WILD
bash ~/Bioinformatics/Selaginella/Shell_scripts/get_counts.sh KRAUS MOEL UNC WILD
GetData -f counts >${Counts}/Selaginella_all.txt
GetData.py -f counts >${Counts}/selag_counts.txt
for species in 'KRAUS' 'MOEL' 'UNC' 'WILD'
do
  GetData.py -f counts -s $species >${Counts}/${species}_counts.txt
done
