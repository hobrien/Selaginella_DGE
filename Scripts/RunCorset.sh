BASEDIR=/Users/HeathOBrien/Bioinformatics

cd $BASEDIR/Mappings
for species in $@
do
  if test $species == 'MOEL'
  then
    corset -g Group1,Group1,Group2,Group2,Group3,Group3 ${species}1_all.bam \
    ${species}1b_all.bam ${species}2_all.bam ${species}2b_all.bam ${species}3_all.bam \
    ${species}3b_all.bam
  else
    corset -g Group1,Group1,Group2,Group2,Group3,Group3,Group4,Group4 ${species}1_all.bam \
    ${species}1b_all.bam ${species}2_all.bam ${species}2b_all.bam ${species}3_all.bam \
    ${species}3b_all.bam ${species}4_all.bam ${species}4b_all.bam
  fi
  mv counts.txt ${species}counts.txt
  mv clusters.txt ${species}clusters.txt
done