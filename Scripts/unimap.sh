for species in $@
do
  echo "testing for ${species} index"
  if ! test -f ~/Index/${species}_nr.1.bt2 || \
     ! test -f ~/Index/${species}_nr.2.bt2 || \
     ! test -f ~/Index/${species}_nr.3.bt2 || \
     ! test -f ~/Index/${species}_nr.4.bt2  || \
     ! test -f ~/Index/${species}_nr.rev.1.bt2 || \
     ! test -f ~/Index/${species}_nr.rev.2.bt2
  then
    bowtie2-build ~/Assemblies/${species}_nr.fa \
    ~/Index/${species}_nr
  fi

  for leaf in $(seq 4)
  do
    if ! test $species == 'MOEL' || ! test $leaf == 4
    then
      echo "mapping ${species}${leaf}"
      if ! test -f ~/Mappings/${species}${leaf}_nr.bam
      then
        bowtie2 -k 1 -p 8 -x ~/Index/${species}_nr \
        -1 <( gunzip -c ~/Reads/${species}${leaf}_1.filtered.fastq.gz ) \
        -2 <( gunzip -c ~/Reads/${species}${leaf}_2.filtered.fastq.gz ) \
        | samtools view -bS - \
        > ~/Mappings/${species}${leaf}_nr.bam
      fi
      if ! test -f ~/Mappings/${species}${leaf}b_nr.bam
      then    
        bowtie2 -k 1 -p 8 -x ~/Index/${species}_nr \
        -1 <( gunzip -c ~/Reads/${species}${leaf}b_1.filtered.fastq.gz ) \
        -2 <( gunzip -c ~/Reads/${species}${leaf}b_2.filtered.fastq.gz ) \
        | samtools view -bS - \
        > ~/Mappings/${species}${leaf}b_nr.bam
      fi
    fi
  done
done
