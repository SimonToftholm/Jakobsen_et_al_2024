## Hisat2 alignment of fastq files

#!bin/bash

for f in $(ls ./*.fastq.gz | rev | cut -c 16- | rev | uniq)
        do
                echo "hisat2 alignment is started for $f"
                hisat2 -p 24 -x /data/home/simontj/Index/Human_Genome/hg38_index --dta -1 ${f}R1_001.fastq.gz -2 ${f}R2_001.fastq.gz -S ${f}.sam
        done
wait


## Samfiles to bamfiles

#!bin/bash

for f in ./*.sam;
        do
                echo "Sam to Bam for $f"
                samtools view -S -b "$f" > "${f%.*}.bam" &
        done
wait

## Sort bamfiles

#!bin/bash

for f in *bam
        do
                echo "Sort $f"
                samtools sort "$f" -o "${f%.*}.sorted.bam" &
        done
wait


## Index bamflies

#!bin/bash

for f in *bam
        do
                echo “Index $f"
                samtools index "$f" &
        done
wait


## Make Count matrix 

featureCounts -p -g gene_id -a gencode.v37.annotation.gtf -T 12 -o count_matrix.txt *sorted.bam 
