BASEDIR=/Users/HeathOBrien/Bioinformatics/Selaginella_DGE

for species in $@
do
  for leaf in $(seq 4)
  do
    if ! test $species == 'MOEL' || ! test $leaf == 4
    then
      adaptor=$(ls $BASEDIR/Raw/Project_Stuart_Casson_S019101b/Sample_${species}${leaf} |head -1 |cut -d_ -f 2)
      fastq-mcf -o <( gunzip -c $BASEDIR/Trimmed/${species}${leaf}_R1.clip.fastq) -o <( gunzip -c $BASEDIR/Trimmed/${species}${leaf}_R2.clip.fastq) $BASEDIR/Adaptors/$adaptor $BASEDIR/Raw/Project_Stuart_Casson_S019101a/Sample_{species}${leaf}/L00*_R1.fq.gz $BASEDIR/Raw/Project_Stuart_Casson_S019101b/Sample_{species}${leaf}/L00*_R2.fq.gz >$BASEDIR/Trimmed/${species}${leaf}.clip.log
      gzip $BASEDIR/Trimmed/${species}${leaf}_R1.clip.fastq $BASEDIR/Trimmed/${species}${leaf}_R2.clip.fastq
      fastq-mcf -o <( gunzip -c $BASEDIR/Trimmed/${species}${leaf}b_R1.clip.fastq) -o <( gunzip -c $BASEDIR/Trimmed/${species}${leaf}b_R2.clip.fastq) $BASEDIR/Adaptors/$adaptor $BASEDIR/Raw/Project_Stuart_Casson_S019101b/Sample_{species}${leaf}/{species}${leaf}_${adaptor}_L00*_R1_001.fastq.gz $BASEDIR/Raw/Project_Stuart_Casson_S019101b/Sample_{species}${leaf}/{species}${leaf}_${adaptor}_L00*_R2_001.fastq.gz >$BASEDIR/Trimmed/${species}${leaf}b.clip.log
      gzip $BASEDIR/Trimmed/${species}${leaf}b_R1.clip.fastq $BASEDIR/Trimmed/${species}${leaf}b_R2.clip.fastq
    fi
  done
done