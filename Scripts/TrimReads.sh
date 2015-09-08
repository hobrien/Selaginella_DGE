BASEDIR=/Users/HeathOBrien/Bioinformatics/Selaginella_DGE

for species in 'KRAUS' 'MOEL' 'UNC' 'WILD'
do
  for leaf in `seq 0 3`
  do
    if  ! test $species = 'MOEL' || ! test $leaf2 = 3 
    then
      gunzip $BASEDIR/Raw/Project_Stuart_Casson_S019101a/Sample_${species}${leaf}/*
      fastq-mcf -o $BASEDIR/Trimmed/${species}${leaf}_R1.clip.fastq \
                -o $BASEDIR/Trimmed/${species}${leaf}_R2.clip.fastq \
                $BASEDIR/Adaptors/${species}${leaf} \
                $BASEDIR/Raw/Project_Stuart_Casson_S019101a/Sample_${species}${leaf}/*_R1.fq \
                $BASEDIR/Raw/Project_Stuart_Casson_S019101a/Sample_${species}${leaf}/*_R2.fq \
                > $BASEDIR/Trimmed/${species}${leaf}.clip.log
      gzip $BASEDIR/Raw/Project_Stuart_Casson_S019101a/Sample_${species}${leaf}/* \
           $BASEDIR/Trimmed/${species}${leaf}_R1.clip.fastq
           $BASEDIR/Trimmed/${species}${leaf}_R2.clip.fastq
      gunzip  $BASEDIR/Raw/Project_Stuart_Casson_S019101b/Sample_${species}${leaf}/*
      fastq-mcf -o $BASEDIR/Trimmed/${species}${leaf}b_R1.clip.fastq \
                -o $BASEDIR/Trimmed/${species}${leaf}b_R2.clip.fastq \
                $BASEDIR/Adaptors/${species}${leaf} \
                $BASEDIR/Raw/Project_Stuart_Casson_S019101b/Sample_${species}${leaf}/*_R1_001.fq \
                $BASEDIR/Raw/Project_Stuart_Casson_S019101b/Sample_${species}${leaf}/*_R2_001.fq \
                > $BASEDIR/Trimmed/${species}${leaf}b.clip.log
      gzip $BASEDIR/Raw/Project_Stuart_Casson_S019101b/Sample_${species}${leaf}/* \
           $BASEDIR/Trimmed/${species}${leaf}b_R1.clip.fastq
           $BASEDIR/Trimmed/${species}${leaf}b_R2.clip.fastq     
    fi
  done
done

