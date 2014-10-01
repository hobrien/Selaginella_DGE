BASEDIR=/Users/HeathOBrien/Bioinformatics
for species in $@
do
  if ! test -f $BASEDIR/Index/${species}_Trinity_all.1.bt2 || \
     ! test -f $BASEDIR/Index/${species}_Trinity_all.2.bt2 || \
     ! test -f $BASEDIR/Index/${species}_Trinity_all.3.bt2 || \
     ! test -f $BASEDIR/Index/${species}_Trinity_all.4.bt2  || \
     ! test -f $BASEDIR/Index/${species}_Trinity_all.rev.1.bt2 || \
     ! test -f $BASEDIR/Index/${species}_Trinity_all.rev.2.bt2
  then
    bowtie2-build $BASEDIR/Assemblies/${species}_Tr.fa \
    $BASEDIR/Index/${species}_Trinity_all
  fi

  for leaf in $(seq 4)
  do
    if ! test $species == 'MOEL' || ! test $leaf == 4
    then
      if ! test -f $BASEDIR/Mappings/${species}${leaf}_all.bam
      then
        bowtie2 -a -x$BASEDIR/Index/${species}_Trinity_all \
        -1 <( gunzip -c $BASEDIR/Reads/${species}${leaf}_1.filtered.fastq.gz ) \
        -2 <( gunzip -c $BASEDIR/Reads/${species}${leaf}_2.filtered.fastq.gz ) \
        | samtools view -bS - \
        > $BASEDIR/Mappings/${species}${leaf}_all.bam
      fi
      if ! test -f $BASEDIR/Mappings/${species}${leaf}b_all.bam
      then    
        bowtie2 -a -x $BASEDIR/Index/${species}_Trinity_all \
        -1 <( gunzip -c $BASEDIR/Reads/${species}${leaf}b_1.filtered.fastq.gz ) \
        -2 <( gunzip -c $BASEDIR/Reads/${species}${leaf}b_2.filtered.fastq.gz ) \
        | samtools view -bS - \
        > $BASEDIR/Mappings/${species}${leaf}b_all.bam
      fi
    fi
  done
done