Index=~/Bioinformatics/Index
Assemblies=~/Bioinformatics/Selaginella/Assemblies
Mappings=~/Bioinformatics/Selaginella/Mappings
Reads=~/Bioinformatics/Selaginella/Reads
GTF=~/Bioinformatics/Selaginella/GTF
Counts=~/Bioinformatics/Selaginella/Counts
for species in $@
do
  if ! test -f ${GTF}/${species}.gtf
  then
    GetData.py -f GTF -s ${species} > ${GTF}/${species}.gtf
  fi
  for leaf in $(seq 4)
  do
    if ! test $species == 'MOEL' || ! test $leaf == 4
    then
      if ! test -f ${Counts}/${species}${leaf}_nr.txt
      then
        if ! test -f ${Mappings}/${species}${leaf}_nr.sam
        then
          samtools view ${Mappings}/${species}${leaf}_nr.bam | sort | perl -p -e 's/_[12]//' \
          > ${Mappings}/${species}${leaf}_nr.sam
        fi
        htseq-count -s no -t CDS ${Mappings}/${species}${leaf}_nr.sam ${GTF}/${species}.gtf \
        > ${Counts}/${species}${leaf}_nr.txt
        #rm ${Mappings}/${species}${leaf}_nr.sam
        AddData.py -f counts -i ${Counts}/${species}${leaf}_nr.txt -n leaf${leaf}
      fi
      if ! test -f ${Counts}/${species}${leaf}b_nr.txt
      then
        if ! test -f ${Mappings}/${species}${leaf}b_nr.sam
        then
          samtools view ${Mappings}/${species}${leaf}b_nr.bam | sort | perl -p -e 's/_[12]//' \
          > ${Mappings}/${species}${leaf}b_nr.sam
        fi
        htseq-count  -s no -t CDS  ${Mappings}/${species}${leaf}b_nr.sam ${GTF}/${species}.gtf \
        > ${Counts}/${species}${leaf}b_nr.txt
        #rm ${Mappings}/${species}${leaf}b_nr.sam
        AddData.py -f counts -i ${Counts}/${species}${leaf}b_nr.txt -n leaf${leaf}b
      fi
    fi
  done
done
