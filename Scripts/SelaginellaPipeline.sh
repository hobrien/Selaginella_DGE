MAXCLUST=30466
BASEDIR=/Users/HeathOBrien/Bioinformatics/Selaginella_DGE

######################################### TRIM AND FILTER READS #####################################
bash $BASEDIR/Scripts/TrimReads.sh

#initialise database with species info
mysql -u root -e "CREATE DATABASE IF NOT EXISTS 'SelaginellaGenomics'"
mysql -u root SelaginellaGenomics < $BASEDIR/Scripts/InitializeSelaginellaGenomics.sql

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

#Add Arabidopsis homologs from ENSEMBL to DB

#Ideally there should be code here to get the homologs. A project for another day
python $BASEDIR/Scripts/AddData.py -f ath -i Homologs/AthHomologs.txt

 ######################################### RUN TRANSDECODER #############################################

#rename amino acid sequences 
for species in 'KRAUS' 'MOEL' 'UNC' 'WILD'
do
  cat $BASEDIR/Transdecoder/${species}_Tr.fa.transdecoder.pep | perl -pe 'BEGIN {$species=shift} s/.*\|g\.(\d+) .*/>${species}|$1/' $species  > $BASEDIR/AA_seqs/${species}.fasta
done

#Add info about coding regions of transcripts to database
for file in `find /Users/HeathOBrien/Bioinformatics/Selaginella/Transdecoder | grep _Tr.fa.transdecoder.bed`
do
	AddData.py -f coding -i $file
done


AddData.py -f ref_coding #reference sequences are all coding, so we don't need info from BED file

#Add orthoMCL info to DB
AddData.py -f orthologs -i /Users/HeathOBrien/Bioinformatics/Selaginella/OrthoMCL/Selaginella_groups.txt

#################################### DO ALL-BY ALL BLAST AND UPLOAD RESULTS ####################################
#Copy reference aa tranlsations sequences AA_seqs
for file in `find $BASEDIR/RefSeq -name *_aa.fa`
do
  if ! test $file == $BASEDIR/RefSeq/Selmo_all_aa.fa
  then
    cat $file | perl -pe 's/scaffold-(\w{4})-(\d+).*/$1|$2/' >$BASEDIR/AA_seqs/$(basename $file | cut -d_ -f1).fasta
  fi
done

cat $BASEDIR/AA_seqs/* >> $BASEDIR/AA_seqs/Selmo_all.fa
makeblastdb -in $BASEDIR/AA_seqs/Selmo_all.fa -dbtype prot -out $BASEDIR/Blast/Selmo_all.fa
blastp -query $BASEDIR/AA_seqs/Selmo_all.fa -db $BASEDIR/Blast/Selmo_all.fa -outfmt 6 -out $BASEDIR/Blast/Selmo_all.fa -num_threads 6

for species in 'KRAUS' 'MOEL' 'UNC' 'WILD'
do
  grep ^$species $BASEDIR/Blast/Selmo_all.bl > $BASEDIR/Blast/$species.bl
  $BASEDIR/Scripts/AddData.py -f blast -i $BASEDIR/Blast/$species.bl
  rm $BASEDIR/Blast/$species.bl
done

#################################### MULTIMAP READS AND RUN CORSET ####################################

#These steps where actually run on BlueCrystal
for species in 'KRAUS' 'MOEL' 'UNC' 'WILD'
do
  bash $BASEDIR/Scripts/multimap.sh
  bash $BASEDIR/Scripts/RunCorset.sh
done

######################## UPLOAD CORSET DATA AND WRITE REPRESENTATIVE SEQS TO FILE ##################
for species in 'KRAUS' 'MOEL' 'UNC' 'WILD'
do
  $BASEDIR/Scripts/AddData.py -f corset_clusters -n $species -i $BASEDIR/Corset/${species}clusters.txt
  $BASEDIR/Scripts/AddData.py -f corset_counts -n $specie -i $BASEDIR/Corset/${species}counts.txt
  $BASEDIR/Scripts/AddData.py -f corset_nr -s $species
  if test $species == 'KRAUS'
  then
    $BASEDIR/Scripts/GetData.py -f corset -s $species -i $BASEDIR/AA_seqs/${species}_aa.fa | perl -pe 's/KRAUS/KRUS/' | perl -pe 's/_/|/' > $BASEDIR/OrthoMCL/compliantFasta/KRUS.fasta
  else
    $BASEDIR/Scripts/GetData.py -f corset -s $species -i $BASEDIR/AA_seqs/${species}_aa.fa | perl -pe 's/_/|/' > $BASEDIR/OrthoMCL/compliantFasta/$species.fasta
  fi
done

################### PLOT LENGTH DIST AND OVERLAP OF NON-REDUNDANT CONTIGS ################
if test -f $BASEDIR/Results/coding_lengths.txt
then
  rm $BASEDIR/Results/coding_lengths.txt
fi

for species in 'KRAUS' 'MOEL' 'UNC' 'WILD'
do
  $BASEDIR/Scripts/GetData.py -f lengths -s $species >> $BASEDIR/Results/coding_lengths.txt
done
echo "coding_lengths.txt: lengths of coding portion of each non-redundant contig" >> $BASEDIR/Results/README.md
Rscript $BASEDIR/Scripts/PlotLengths.R $BASEDIR/Results/coding_lengths.txt $BASEDIR/Figures/coding_lengths.png
echo "coding_lengths.png: histograms of lengths of coding portions for non-redundant contig" >> $BASEDIR/Figures/README.md

###### IDENTIFY NON-REDUNDANT CONTIGS WITH > 1 CODING SEQUENCE AND BLAST AGAINST REF #####

if ! test -d $BASEDIR/Chimeras
then
  mkdir $BASEDIR/Chimeras
fi
if ! test -f $BASEDIR/Blast/Selmo_ensembl.psq || if ! test -f $BASEDIR/Blast/Selmo_ensembl.pin || if ! test -f $BASEDIR/Blast/Selmo_ensembl.phr
then
  makeblastdb -in $BASEDIR/RefSeq/Selaginella_moellendorffii.v1.0.17.pep.all.fa -dbtype prot -out $BASEDIR/Blast/Selmo_ensembl
fi

for species in 'KRAUS' 'MOEL' 'UNC' 'WILD'
do
  $BASEDIR/Scripts/GetData.py -f chimeric -s $species > $BASEDIR/Chimeras/$species.fa
  blastx -query $BASEDIR/Chimeras/$species.fa -db $BASEDIR/Blast/Selmo_ensembl -outfmt 6 -out $BASEDIR/Chimeras/$species.bl -num_threads 6
  $BASEDIR/Scripts/AddData.py -f blastx -i $BASEDIR/Chimeras/$species.bl
done

######################################### RUN ORTHOMCL #############################################
cd $BASEDIR/OrthoMCL/
orthomclFilterFasta $BASEDIR/OrthoMCL/compliantFasta/ 10 20 #I'm not sure if this file is actually needed for anything
$BASEDIR/Scripts/FilterBlast.py -i $BASEDIR/Blast/Selmo_all.bl > $BASEDIR/OrthoMCL/goodProteins.bl
mysql -u root -e "CREATE DATABASE IF NOT EXISTS 'orthomcl'"
mysql -u root orthomcl < $BASEDIR/Scripts/InitializeDB.sql
orthomclInstallSchema $BASEDIR/OrthoMCL/orthomcl.config.template
orthomclBlastParser $BASEDIR/OrthoMCL/goodProteins.bl $BASEDIR/OrthoMCL/compliantFasta > $BASEDIR/OrthoMCL/similarSequences.txt
orthomclLoadBlast $BASEDIR/OrthoMCL/orthomcl.config.template $BASEDIR/OrthoMCL/similarSequences.txt
orthomclPairs $BASEDIR/OrthoMCL/orthomcl.config.template $BASEDIR/OrthoMCL/orthomclPair.log cleanup=yes
rm -r $BASEDIR/OrthoMCL/pairs
orthomclDumpPairsFiles orthomcl.config.template
mcl $BASEDIR/OrthoMCL/mclInput --abc -I 1.2 -o $BASEDIR/OrthoMCL/MCL/out.data.mci.I12
orthomclMclToGroups OG2_ 0 < $BASEDIR/OrthoMCL/MCL/out.data.mci.I12 >$BASEDIR/OrthoMCL/Selaginella_groups.txt
mysql -u root SelaginellaGenomics < $BASEDIR/Scripts/DeleteOrthologs.sql
$BASEDIR/Scripts/AddData.py -f orthologs -i $BASEDIR/OrthoMCL/Selaginella_groups.txt

#Add info about which gene families are plastid encoded
echo "UPDATE OrthoGroups, OrthoGroupInfo SET OrthoGroupInfo.genome = 'plastid' WHERE OrthoGroups.orthoID = OrthoGroupInfo.orthoID AND OrthoGroups.geneID LIKE 'EFJ_ADH%';" | mysql -u root SelaginellaGenomics

############################# MAKE VENN DIAGRAM OF ORTHOLOG OVERLAPS ###############################
for species in 'KRAUS' 'MOEL' 'UNC' 'WILD'
do
  $BASEDIR/Scripts/GetData.py -f clusters -s $species > $BASEDIR/Results/OrthoGroups/$species.txt
done
echo "OrthoGroups: lists of ortholog groups for each species" >> $BASEDIR/Results/README.md
cd $BASEDIR/Figures
$BASEDIR/Scripts/CompareSets.pl $BASEDIR/Results/OrthoGroups/UNC.txt \
                                $BASEDIR/Results/OrthoGroups/MOEL.txt \
                                $BASEDIR/Results/OrthoGroups/KRAUS.txt \
                                $BASEDIR/Results/OrthoGroups/WILD.txt 
mv VennPlot.png $BASEDIR/Figures/ortholog_overlap.png
echo "ortholog_overlap.png: venn diagram of ortholog group overlaps" >> $BASEDIR/Figures/README.md

############################## DIFFERENTIAL EXPRESSION ANALYSES ###########################
for species in 'KRAUS' 'MOEL' 'UNC' 'WILD'
do
  $BASEDIR/Scripts/GetData.py -f count_totals -s $species > $BASEDIR/Corset/${species}count_totals.txt
  #$BASEDIR/Scripts/FilterCounts.py --in $BASEDIR/Corset/${species}count_totals.txt --min 30 >$BASEDIR/DGEclust/${species}counts.txt
  $BASEDIR/Scripts/RunDGEclust.py  $BASEDIR/Corset/${species}count_totals.txt
  for leaf1 in `seq 0 2`
  do
    for leaf2 in `seq $(($leaf1+1)) 3`
    do
      if ! test $species = 'MOEL' || ! test $leaf2 = 3
      then
        python $BASEDIR/Scripts/GetDGEclustPvals.py $BASEDIR/Corset/KRAUScount_totals.txt leaf1 leaf2
        Rscript $BASEDIR/Scripts/DESeq2_v_DGEClust.R $BASEDIR/DGEclust/${species}_$(($leaf1+1))$(($leaf2+1)).txt $BASEDIR/DGEclust/${species}counts.txt 0.01
        python $BASEDIR/Scripts/AddData.py -f de -n ${species} -i $BASEDIR/DGEclust/${species}_$(($leaf1+1))$(($leaf2+1))_0.01_overlap2.txt
      fi
    done
  done
done
Rscript $BASEDIR/Scripts/VolcanoPlot.R $BASEDIR/Figures/Volcano.png

############################# MAKE LIST OF CHIMERIC DE GENES ###############################

rm /tmp/chimeras.txt
mysql -u root SelaginellaGenomics < $BASEDIR/Scripts/DEchimeras.sql
cp /tmp/chimeras.txt $BASEDIR/Results/chimeras.txt
echo "chimeras.txt: list of DE genes with 2 or more CDSs" >> $BASEDIR/Results/README.md

############################# MAKE VENN DIAGRAM OF DE GENES ###############################
if ! test -d $BASEDIR/Results/DEgroups/
then
  mkdir $BASEDIR/Results/DEgroups/
fi
for species in 'KRAUS' 'MOEL' 'UNC' 'WILD'
do
  $BASEDIR/Scripts/GetData.py -f de -s $species > $BASEDIR/Results/DEgroups/$species.txt
done
echo "DEgroups: lists of ortholog groups with at least 1 member DE in at least 1 comparison for each species" >> $BASEDIR/Results/README.md
cd $BASEDIR/Figures
$BASEDIR/Scripts/CompareSets.pl $BASEDIR/Results/DEgroups/UNC.txt \
                                $BASEDIR/Results/DEgroups/MOEL.txt \
                                $BASEDIR/Results/DEgroups/KRAUS.txt \
                                $BASEDIR/Results/DEgroups/WILD.txt 
mv VennPlot.png $BASEDIR/Figures/de_ortholog_overlap.png
echo "de_ortholog_overlap.png: venn diagram of overlap of DE genes" >> $BASEDIR/Figures/README.md

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

################################ MAP READS TO NR SEQS AND GET HIT COUNTS ##################################
cd /Users/HeathOBrien/Bioinformatics/Selaginella/Non_redundant
for species in 'KRAUS' 'MOEL' 'UNC' 'WILD'
do
  GetData.py -f nr_seqs -s $species > ${species}_nr.fa
  GetData.py -f nr_seq -s $species > ${species}
done
CompareSets.pl UNC MOEL KRAUS WILD
mv VennPlot.pdf ~/Google\ Drive/Selaginella/Figures/OrthoVenn.pdf
  
bash ~/Bioinformatics/Selaginella/Shell_scripts/unimap.sh KRAUS MOEL UNC WILD #this step was run on BlueCrystal
bash ~/Bioinformatics/Selaginella/Shell_scripts/get_counts.sh KRAUS MOEL UNC WILD
GetData -f counts >${Counts}/Selaginella_all.txt
GetData.py -f counts >${Counts}/selag_counts.txt
for species in 'KRAUS' 'MOEL' 'UNC' 'WILD'
do
  GetData.py -f counts -s $species >${Counts}/${species}_counts.txt
done
