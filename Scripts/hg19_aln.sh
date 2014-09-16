bowtie2 -p 8 -x hg19 -1 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/UNC1b_R1.clip.fastq -2 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/UNC1b_R2.clip.fastq -S /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC1b_hg19.sam
perl -pi -e 's/ (\d).*/_$1/' /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC1b_hg19.sam
samtools view -bS /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC1b_hg19.sam > /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC1b_hg19.bam
bam2fastq --no-aligned -o /Users/HeathOBrien/Bioinformatics/Selaginella/Filtered/UNC1b\#.filtered.fastq /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC1b_hg19.bam 

bowtie2 -p 8 -x hg19 -1 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/UNC2b_R1.clip.fastq -2 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/UNC2b_R2.clip.fastq -S /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC2b_hg19.sam
perl -pi -e 's/ (\d).*/_$1/' /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC2b_hg19.sam
samtools view -bS /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC2b_hg19.sam >/Users/HeathOBrien/Bioinformatics/Selaginella/ >Mappings/UNC2b_hg19.bam
bam2fastq --no-aligned -o /Users/HeathOBrien/Bioinformatics/Selaginella/Filtered/UNC2b\#.filtered.fastq /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC2b_hg19.bam 

bowtie2 -p 8 -x hg19 -1 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/UNC3b_R1.clip.fastq -2 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/UNC3b_R2.clip.fastq -S /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC3b_hg19.sam
perl -pi -e 's/ (\d).*/_$1/' /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC3b_hg19.sam
samtools view -bS /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC3b_hg19.sam >/Users/HeathOBrien/Bioinformatics/Selaginella/ >Mappings/UNC3b_hg19.bam
bam2fastq --no-aligned -o /Users/HeathOBrien/Bioinformatics/Selaginella/Filtered/UNC3b\#.filtered.fastq /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC3b_hg19.bam 

bowtie2 -p 8 -x hg19 -1 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/UNC4b_R1.clip.fastq -2 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/UNC4b_R2.clip.fastq -S /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC4b_hg19.sam
perl -pi -e 's/ (\d).*/_$1/' /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC4b_hg19.sam
samtools view -bS /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC4b_hg19.sam >/Users/HeathOBrien/Bioinformatics/Selaginella/ >Mappings/UNC4b_hg19.bam
bam2fastq --no-aligned -o /Users/HeathOBrien/Bioinformatics/Selaginella/Filtered/UNC4b\#.filtered.fastq /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/UNC4b_hg19.bam 

bowtie2 -p 8 -x hg19 -1 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/WILD1b_R1.clip.fastq -2 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/WILD1b_R2.clip.fastq -S /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/WILD1b_hg19.sam
perl -pi -e 's/ (\d).*/_$1/' /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/WILD1b_hg19.sam
samtools view -bS /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/WILD1b_hg19.sam >/Users/HeathOBrien/Bioinformatics/Selaginella/ >Mappings/WILD1b_hg19.bam
bam2fastq --no-aligned -o /Users/HeathOBrien/Bioinformatics/Selaginella/Filtered/WILD1b\#.filtered.fastq /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/WILD1b_hg19.bam 

bowtie2 -p 8 -x hg19 -1 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/WILD2b_R1.clip.fastq -2 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/WILD2b_R2.clip.fastq -S /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/WILD2b_hg19.sam
perl -pi -e 's/ (\d).*/_$1/' /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/WILD2b_hg19.sam
samtools view -bS /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/WILD2b_hg19.sam >/Users/HeathOBrien/Bioinformatics/Selaginella/ >Mappings/WILD2b_hg19.bam
bam2fastq --no-aligned -o /Users/HeathOBrien/Bioinformatics/Selaginella/Filtered/WILD2b\#.filtered.fastq /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/WILD2b_hg19.bam 

bowtie2 -p 8 -x hg19 -1 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/WILD3b_R1.clip.fastq -2 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/WILD3b_R2.clip.fastq -S /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/WILD3b_hg19.sam
perl -pi -e 's/ (\d).*/_$1/' /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/WILD3b_hg19.sam
samtools view -bS /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/WILD3b_hg19.sam >/Users/HeathOBrien/Bioinformatics/Selaginella/ >Mappings/WILD3b_hg19.bam
bam2fastq --no-aligned -o /Users/HeathOBrien/Bioinformatics/Selaginella/Filtered/WILD3b\#.filtered.fastq /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/WILD3b_hg19.bam 

bowtie2 -p 8 -x hg19 -1 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/WILD4b_R1.clip.fastq -2 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/WILD4b_R2.clip.fastq -S /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/WILD4b_hg19.sam
perl -pi -e 's/ (\d).*/_$1/' /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/WILD4b_hg19.sam
samtools view -bS /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/WILD4b_hg19.sam >/Users/HeathOBrien/Bioinformatics/Selaginella/ >Mappings/WILD4b_hg19.bam
bam2fastq --no-aligned -o /Users/HeathOBrien/Bioinformatics/Selaginella/Filtered/WILD4b\#.filtered.fastq /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/WILD4b_hg19.bam 

bowtie2 -p 8 -x hg19 -1 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/KRAUS1b_R1.clip.fastq -2 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/KRAUS1b_R2.clip.fastq -S /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/KRAUS1b_hg19.sam
perl -pi -e 's/ (\d).*/_$1/' /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/KRAUS1b_hg19.sam
samtools view -bS /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/KRAUS1b_hg19.sam >/Users/HeathOBrien/Bioinformatics/Selaginella/ >Mappings/KRAUS1b_hg19.bam
bam2fastq --no-aligned -o /Users/HeathOBrien/Bioinformatics/Selaginella/Filtered/KRAUS1b\#.filtered.fastq /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/KRAUS1b_hg19.bam 

bowtie2 -p 8 -x hg19 -1 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/KRAUS2b_R1.clip.fastq -2 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/KRAUS2b_R2.clip.fastq -S /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/KRAUS2b_hg19.sam
perl -pi -e 's/ (\d).*/_$1/' /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/KRAUS2b_hg19.sam
samtools view -bS /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/KRAUS2b_hg19.sam >/Users/HeathOBrien/Bioinformatics/Selaginella/ >Mappings/KRAUS2b_hg19.bam
bam2fastq --no-aligned -o /Users/HeathOBrien/Bioinformatics/Selaginella/Filtered/KRAUS2b\#.filtered.fastq /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/KRAUS2b_hg19.bam 

bowtie2 -p 8 -x hg19 -1 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/KRAUS3b_R1.clip.fastq -2 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/KRAUS3b_R2.clip.fastq -S /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/KRAUS3b_hg19.sam
perl -pi -e 's/ (\d).*/_$1/' /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/KRAUS3b_hg19.sam
samtools view -bS /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/KRAUS3b_hg19.sam >/Users/HeathOBrien/Bioinformatics/Selaginella/ >Mappings/KRAUS3b_hg19.bam
bam2fastq --no-aligned -o /Users/HeathOBrien/Bioinformatics/Selaginella/Filtered/KRAUS3b\#.filtered.fastq /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/KRAUS3b_hg19.bam 

bowtie2 -p 8 -x hg19 -1 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/KRAUS4b_R1.clip.fastq -2 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/KRAUS4b_R2.clip.fastq -S /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/KRAUS4b_hg19.sam
perl -pi -e 's/ (\d).*/_$1/' /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/KRAUS4b_hg19.sam
samtools view -bS /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/KRAUS4b_hg19.sam >/Users/HeathOBrien/Bioinformatics/Selaginella/ >Mappings/KRAUS4b_hg19.bam
bam2fastq --no-aligned -o /Users/HeathOBrien/Bioinformatics/Selaginella/Filtered/KRAUS4b\#.filtered.fastq /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/KRAUS4b_hg19.bam 

bowtie2 -p 8 -x hg19 -1 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/MOEL1b_R1.clip.fastq -2 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/MOEL1b_R2.clip.fastq -S /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/MOEL1b_hg19.sam
perl -pi -e 's/ (\d).*/_$1/' /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/MOEL1b_hg19.sam
samtools view -bS /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/MOEL1b_hg19.sam >/Users/HeathOBrien/Bioinformatics/Selaginella/ >Mappings/MOEL1b_hg19.bam
bam2fastq --no-aligned -o /Users/HeathOBrien/Bioinformatics/Selaginella/Filtered/MOEL1b\#.filtered.fastq /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/MOEL1b_hg19.bam 

bowtie2 -p 8 -x hg19 -1 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/MOEL2b_R1.clip.fastq -2 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/MOEL2b_R2.clip.fastq -S /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/MOEL2b_hg19.sam
perl -pi -e 's/ (\d).*/_$1/' /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/MOEL2b_hg19.sam
samtools view -bS /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/MOEL2b_hg19.sam >/Users/HeathOBrien/Bioinformatics/Selaginella/ >Mappings/MOEL2b_hg19.bam
bam2fastq --no-aligned -o /Users/HeathOBrien/Bioinformatics/Selaginella/Filtered/MOEL2b\#.filtered.fastq /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/MOEL2b_hg19.bam 

bowtie2 -p 8 -x hg19 -1 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/MOEL3b_R1.clip.fastq -2 /Users/HeathOBrien/Bioinformatics/Selaginella/Trimmed/MOEL3b_R2.clip.fastq -S /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/MOEL3b_hg19.sam
perl -pi -e 's/ (\d).*/_$1/' /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/MOEL3b_hg19.sam
samtools view -bS /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/MOEL3b_hg19.sam >/Users/HeathOBrien/Bioinformatics/Selaginella/ >Mappings/MOEL3b_hg19.bam
bam2fastq --no-aligned -o /Users/HeathOBrien/Bioinformatics/Selaginella/Filtered/MOEL3b\#.filtered.fastq /Users/HeathOBrien/Bioinformatics/Selaginella/Mappings/MOEL3b_hg19.bam 
