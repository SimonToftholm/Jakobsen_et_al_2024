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


## Call peaks with MACS2 

#!bin/bash

for f in *sorted.dup.bam;
        do
                echo "MACS2 peak calling for $f"
                macs2 callpeak -t "$f"-f BAM -g hs -n ${f%.*} -q 0.05
        done


## Index barflies

#!bin/bash

for f in *bam
        do
                echo “Index $f"
                samtools index "$f" &
        done
wait


## Make BigWig files

#!bin/bash

for f in *.bam
        do
                bamCoverage -b "$f" --normalizeUsing RPKM -o "${f%.*}.bw" -p 4 &
        done
wait


## Make TagDirectories

#!bin/bash 

for SAMFILE in *.bam

        do
                echo "Make Tagdirectory on " $SAMFILE
                # Make tagdirectory and capture output
                 makeTagDirectory $SAMFILE.TD $SAMFILE -keepAll -sspe -genome hg38 -checkGC -tbp 1 2 -single -fragLength 200 &
        done
