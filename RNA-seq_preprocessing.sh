{\rtf1\ansi\ansicpg1252\cocoartf2709
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;\f1\fnil\fcharset0 Menlo-Regular;\f2\fnil\fcharset0 Menlo-Bold;
}
{\colortbl;\red255\green255\blue255;\red56\green185\blue199;\red0\green0\blue0;\red57\green192\blue38;
\red86\green32\blue244;\red219\green39\blue218;\red170\green171\blue37;\red202\green51\blue35;}
{\*\expandedcolortbl;;\cssrgb\c25546\c77007\c82023;\csgray\c0;\cssrgb\c25706\c77963\c19557;
\cssrgb\c41681\c25958\c96648;\cssrgb\c89513\c29736\c88485;\cssrgb\c72331\c71682\c18599;\cssrgb\c83899\c28663\c18026;}
\paperw11900\paperh16840\margl1440\margr1440\vieww51940\viewh31760\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs32 \cf0 ## Hisat2 alignment of fastq files\
\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f1\fs26 \cf2 \CocoaLigature0 #!bin/bash\cf3 \
\
\cf4 for\cf3  f \cf4 in\cf3  \cf4 $(\cf3 ls ./*.fastq.gz \cf4 |\cf3  rev \cf4 |\cf3  
\f2\b \cf5 cut\cf6  -c
\f1\b0 \cf3  16- \cf4 |\cf3  rev \cf4 |\cf3  uniq\cf4 )\cf3 \
        \cf4 do\cf3 \
                
\f2\b \cf5 echo
\f1\b0 \cf3  
\f2\b \cf7 "hisat2 alignment is started for $f"
\f1\b0 \cf3 \
                hisat2
\f2\b \cf6  -p
\f1\b0 \cf3  24
\f2\b \cf6  -x
\f1\b0 \cf3  /data/home/simontj/Index/Human_Genome/hg38_index
\f2\b \cf6  --dta
\f1\b0 \cf3  -1 
\f2\b \cf8 $\{f\}
\f1\b0 \cf3 R1_001.fastq.gz -2 
\f2\b \cf8 $\{f\}
\f1\b0 \cf3 R2_001.fastq.gz
\f2\b \cf6  -S
\f1\b0 \cf3  
\f2\b \cf8 $\{f\}
\f1\b0 \cf3 .sam\
        \cf4 done\cf3 \
wait\
\
\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs32 \cf0 \CocoaLigature1 ## Samfiles to bamfiles\
\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f1\fs26 \cf2 \CocoaLigature0 #!bin/bash\cf3 \
\
\cf4 for\cf3  f \cf4 in\cf3  ./*.sam\cf4 ;\cf3 \
        \cf4 do\cf3 \
                
\f2\b \cf5 echo
\f1\b0 \cf3  
\f2\b \cf7 "Sam to Bam for $f"
\f1\b0 \cf3 \
                samtools view
\f2\b \cf6  -S -b
\f1\b0 \cf3  
\f2\b \cf7 "$f"
\f1\b0 \cf3  \cf4 >\cf3  
\f2\b \cf7 "$\{f%.*\}.bam"
\f1\b0 \cf3  \cf4 &\cf3 \
        \cf4 done\cf3 \
wait\
\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs32 \cf0 \CocoaLigature1 ## Sort bamfiles\
\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f1\fs26 \cf2 \CocoaLigature0 #!bin/bash\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f0\fs32 \cf0 \CocoaLigature1 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f1\fs26 \cf4 \CocoaLigature0 for\cf3  f \cf4 in\cf3  *bam\
        \cf4 do\cf3 \
                
\f2\b \cf5 echo
\f1\b0 \cf3  
\f2\b \cf7 "Sort $f"
\f1\b0 \cf3 \
                samtools 
\f2\b \cf5 sort
\f1\b0 \cf3  
\f2\b \cf7 "$f"\cf6  -o
\f1\b0 \cf3  
\f2\b \cf7 "$\{f%.*\}.sorted.bam"
\f1\b0 \cf3  \cf4 &\cf3 \
        \cf4 done\cf3 \
wait
\f0\fs32 \cf0 \CocoaLigature1 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f1\fs26 \cf3 \CocoaLigature0 \
\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs32 \cf0 \CocoaLigature1 ## Index bamflies\
\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f1\fs26 \cf2 \CocoaLigature0 #!bin/bash\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f0\fs32 \cf0 \CocoaLigature1 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f1\fs26 \cf4 \CocoaLigature0 for\cf3  f \cf4 in\cf3  *bam\
        \cf4 do\cf3 \
                
\f2\b \cf5 echo
\f1\b0 \cf3  
\f2\b \cf7 \'93Index $f"
\f1\b0 \cf3 \
                samtools 
\f2\b \cf5 index
\f1\b0 \cf3  
\f2\b \cf7 "$f"
\f1\b0 \cf3  \cf4 &\cf3 \
        \cf4 done\cf3 \
wait\
\
\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs32 \cf0 \CocoaLigature1 ## Make Count matrix 
\f1\fs26 \cf3 \CocoaLigature0 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\
featureCounts -p -g gene_id -a gencode.v37.annotation.gtf -T 12 -o count_matrix.txt *sorted.bam }