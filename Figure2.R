##############################################
################## FIGURE 2 ##################
##############################################

################# Packages ################# 
library(tidyverse)
library(mltools)
library(fgsea)
library(msigdbr)
library(gridExtra)
library(cBioPortalData)
library(ggplot2)
library(survminer)
library(survival)
library(ggsurvfit)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(plyr)

################# Figure 2A ################## 
## Count reads at MYC concensus regions in BT549 using homer on cluster
# annotatePeaks.pl ../MYC_Project/MYC_MED1_ChIP/peak/MYC_BT549_Consensus.bed hg38 -d ../MYC_Project/MYC_MED1_ChIP/tagdir/Combined_MYC_BT549/ -size 1000)
BT549_count <- read.delim("MYC_BT549_Consensus.count.txt")

# Identify MYC promoters and enhancers
BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]

# Summarise enhancer contribution 
BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Gene.Name","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_BT549..Tag.Count.in.1000.bp..40939315.0.Total..normalization.factor...0.24..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

# Summarise promoter contribution 
BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Gene.Name","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_BT549..Tag.Count.in.1000.bp..40939315.0.Total..normalization.factor...0.24..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

# Merge genes with either promoter or enhancer MYC contribution
plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0

# Calculate the difference en contribution from different regions
plot$count <- plot$Count_prom-plot$Count_enh

# Add gene expression information to the genes with nearby MYC binding if the gene is significant regulated upon siRNA knockdown (data from DESeq2 script)
plot2 <- read.table("DE_analysis_KO_BT549.txt")
plot2 <- plot2[plot2$padj_MYC_vs_Crtl < .05,]
plot2 <- plot2[,c("symbol", "Log2FC_MYC_vs_Crtl")]
plot3 <- left_join(plot, plot2)
plot3[is.na(plot3)] = 0
plot3 <- plot3[!duplicated(plot3$symbol),]

# Add HI-C data from ABC-analysis
# ABC analysis

# The ABC analysis was run on the cluster following the instruction from the paper with HIC data 

## Call candidate enhancers (public ATAC-seq and H3K27ac ChIP-seq from MCF7 and BT549 cells were used including gene expression data from DepMap)
# 1. python ABC-Enhancer-Gene-Prediction/src/makeCandidateRegions.py --narrowPeak BT549_peaks.narrowPeak.sorted --bam BT549_ATAC_SRR13755451_SRR13755452_sorted.bam --outDir BT549_Candidate_Peak --chrom_sizes hg38.chromsizes.new --regions_blocklist Blacklist_hg38.bed --regions_includelist ABC-Enhancer-Gene-Prediction/reference/hg38/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed --peakExtendFromSummit 250 --nStrongestPeaks 15000

## Determine enhancer activity
# 2. python ABC-Enhancer-Gene-Prediction/src/run.neighborhoods.py --candidate_enhancer_regions BT549_Candidate_Peak/BT549_peaks.narrowPeak.sorted.candidateRegions.bed --genes ABC-Enhancer-Gene-Prediction/reference/hg38/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.bed --H3K27ac BT549_H3K27ac_sorted.bam --DHS BT549_ATAC_SRR13755451_SRR13755452_sorted.bam --expression_table BT549_counts_final.txt --chrom_sizes hg38.chromsizes.new --outdir BT549_Neighborhoods/

## Prepare HIC data with Juicer (public HI-C data from MCF7, BT549 and TNBC patients)
# 3. /data/home/simontj/Tools/juicer/CPU/juicer.sh -z /data/home/simontj/Index/reference/hg38.fa -p /data/home/simontj/Rasmus/MYC_Project/ABC/hg38.chromsizes -D /data/home/simontj/Tools/juicer/CPU/aidenlab/ -g hg38 -t 40 &
# 4. ABC-Enhancer-Gene-Prediction/src/juicebox_dump.py --hic_file HIC_public/BT549/aligned/inter_30.hic --juicebox "java -jar Juicebox_1.11.08.jar" --outdir BT549_HIC &

## Calculate ABC score 
# 5. python ABC-Enhancer-Gene-Prediction/src/predict.py --enhancers BT549_Neighborhoods/EnhancerList.txt --genes BT549_Neighborhoods/GeneList.txt --HiCdir BT549_HIC_Res5000/ --chrom_sizes hg38.chromsizes.new --hic_resolution 5000 --scale_hic_using_powerlaw --threshold .05 --outdir BT549_Predictions/

# Add ABC score to the matrix for ranking
ABC_BT549 <- read.delim("BT549_EnhancerPredictionsFull.txt")

# Run on cluster to count MYC binding at ABC used enhancer for subsetting
## annotatePeaks.pl EnhancerPredictionsFull.txt hg38 -d ../../CellLine_MYC_MED1_ChIP/tagdir/Combined_MYC_BT549/ > EnhancerPredictionsFull_MYC_binding.txt
ABC_BT549_MYC <- read.delim("EnhancerPredictionsFull_MYC_binding.txt")

# Add information to ABC dataframe
ABC_merged <- merge(ABC_BT549, ABC_BT549_MYC, by.x = "name", by.y = "PeakID..cmd.annotatePeaks.pl.EnhancerPredictionsFull.txt.hg38..d.......CellLine_MYC_MED1_ChIP.tagdir.Combined_MYC_BT549..")

# Subset connection with MYC binding
ABC_merged <- ABC_merged[ABC_merged$......CellLine_MYC_MED1_ChIP.tagdir.Combined_MYC_BT549..Tag.Count.in.given.bp..40939315.0.Total..normalization.factor...0.24..effective.total...10000000. > 30,]

ABC_plot <- ABC_merged[,c("TargetGene","ABC.Score","powerlaw_contact")]
ABC_plot %>% group_by(TargetGene) %>% summarise_all(funs(sum)) -> ABC_plot_1
colnames(ABC_plot_1) <- c("symbol","ABC.Score","HIC")


plot4 <- left_join(plot3, ABC_plot_1)
plot4[is.na(plot4)] = 0

# Summarize bins for plotting
plot5 <- plot4[1:12450,2:7]
plot5 <- plot5[order(plot5$count),]
plot6 <- rowsum(plot5,rep(1:249,each=50))

ggplot() +
  geom_smooth(data= plot6 ,mapping=aes(x=rank(-count),y=Count_enh/50)) +
  geom_smooth(data= plot6 ,mapping=aes(x=rank(-count),y=Count_prom/50)) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + xlab("Ranked Promoter - Enhancer dominace") + ylab("MYC tag density") +
  theme(plot.title = element_text(hjust = 0.5))



# Density for ABC sscore over enhancer promoter ranking
ggplot() +
  geom_smooth(data= plot6 ,mapping=aes(x=rank(-count),y=(ABC.Score/50)), size = 2) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + xlab("Ranked genes") + ylab("ABC score") +
  theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(ylim=c(0,0.17))



################# Figure 2B ################## 
# Bigwig files were loaded into IGV genome browser software and the region of interest were captured

################# Figure 2C ################## 
# Gene expression over the time frame
# Subset genes which are upregulated by MYC (significant genes are only including from earlier creation of the matrix padj<0.05) 
Up <- plot4[plot4$Log2FC_MYC_vs_Crtl > 0, ]

# Bin genes based on enhancer and promoter dependency
Up[, "group"] <- bin_data(Up$count, bins=4, binType = "quantile")
one.way <- aov(abs(Log2FC_MYC_vs_Crtl) ~ group, data = Up)
summary(one.way)

# Plot bins
ggplot() +
  geom_boxplot(data=Up, mapping=aes(x=fct_rev(group), y=(Log2FC_MYC_vs_Crtl)), alpha=0)+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + xlab("Ranked Promoter - Enhancer dominace") + ylab("Gene Log2fc(siMYC/siCrlt)") +
  theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(ylim=c(0,3)) 

Down <- plot4[plot4$Log2FC_MYC_vs_Crtl < 0,]

# Subset genes which are downregulated by MYC (significant genes are only including from earlier creation of the matrix padj<0.05)
Down[, "group"] <- bin_data(Down$count, bins=4, binType = "quantile")

# Bin genes based on enhancer and promoter dependency
one.way <- aov(abs(Log2FC_MYC_vs_Crtl) ~ group, data = Down)
summary(one.way)

# Plot bins
ggplot() +
  geom_boxplot(data=Down, mapping=aes(x=fct_rev(group), y=(Log2FC_MYC_vs_Crtl)), alpha=0) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + xlab("Ranked Promoter - Enhancer dominace") + ylab("Gene Log2fc(siMYC/siCrlt)") +
  theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(ylim=c(-1.5,0))


################# Figure 2D ################## 
# Load in gene sets for enrichment analysis
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

# MCF7 enhancer promoter ranking
BT549_count <- read.delim("MYC_MCF7_Consensus.count.txt")
BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_MC7..Tag.Count.in.1000.bp..42607752.5.Total..normalization.factor...0.23..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_MC7..Tag.Count.in.1000.bp..42607752.5.Total..normalization.factor...0.23..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

MCF7_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
MCF7_prom <- as.data.frame(MCF7_prom@result)

MCF7_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
MCF7_enh <- as.data.frame(MCF7_enh@result)



# BT549 enhancer promoter ranking
BT549_count <- read.delim("MYC_BT549_Consensus.count.txt")
BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_BT549..Tag.Count.in.1000.bp..40939315.0.Total..normalization.factor...0.24..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_BT549..Tag.Count.in.1000.bp..40939315.0.Total..normalization.factor...0.24..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

BT549_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
BT549_prom <- as.data.frame(BT549_prom@result)

BT549_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g) 
BT549_enh <- as.data.frame(BT549_enh@result)

# H2171 enhancer promoter ranking
BT549_count <- read.delim("MYC_H2171_enhancer.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","MYC_H2171..Tag.Count.in.1000.bp..30710431.0.Total..normalization.factor...0.33..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","MYC_H2171..Tag.Count.in.1000.bp..30710431.0.Total..normalization.factor...0.33..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

H2171_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
H2171_prom <- as.data.frame(H2171_prom@result)

H2171_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
H2171_enh <- as.data.frame(H2171_enh@result)





# U2OS enhancer promoter ranking
BT549_count <- read.delim("MYC_U2OS_enhancer.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","MYC_U2OS.Tag.Count.in.1000.bp..5881570.0.Total..normalization.factor...1.70..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","MYC_U2OS.Tag.Count.in.1000.bp..5881570.0.Total..normalization.factor...1.70..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

U2OS_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
U2OS_prom <- as.data.frame(U2OS_prom@result)

U2OS_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
U2OS_enh <- as.data.frame(U2OS_enh@result)


# H128 enhancer promoter ranking
BT549_count <- read.delim("MYC_H128_enhancer.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","MYC_H128.Tag.Count.in.1000.bp..16989736.0.Total..normalization.factor...0.59..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","MYC_H128.Tag.Count.in.1000.bp..16989736.0.Total..normalization.factor...0.59..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

H128_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
H128_prom <- as.data.frame(H128_prom@result)

H128_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
H128_enh <- as.data.frame(H128_enh@result)

# HT1080 enhancer promoter ranking
BT549_count <- read.delim("MYC_HT1080_enhancer.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","MYC_HT1080.Tag.Count.in.1000.bp..15495594.0.Total..normalization.factor...0.65..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","MYC_HT1080.Tag.Count.in.1000.bp..15495594.0.Total..normalization.factor...0.65..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

HT1080_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
HT1080_prom <- as.data.frame(HT1080_prom@result)

HT1080_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
HT1080_enh <- as.data.frame(HT1080_enh@result)

# Hela enhancer promoter ranking
BT549_count <- read.delim("MYC_HeLA_enhancer.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","MYC_HeLA.Tag.Count.in.1000.bp..12570053.0.Total..normalization.factor...0.80..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","MYC_HeLA.Tag.Count.in.1000.bp..12570053.0.Total..normalization.factor...0.80..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

HeLa_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
HeLa_prom <- as.data.frame(HeLa_prom@result)

HeLa_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
HeLa_enh <- as.data.frame(HeLa_enh@result)

# MM1S enhancer promoter ranking
BT549_count <- read.delim("MYC_MM1S_enhancer.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","MYC_MM1S.Tag.Count.in.1000.bp..19719018.0.Total..normalization.factor...0.51..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","MYC_MM1S.Tag.Count.in.1000.bp..19719018.0.Total..normalization.factor...0.51..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

MM1S_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
MM1S_prom <- as.data.frame(MM1S_prom@result)

MM1S_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
MM1S_enh <- as.data.frame(MM1S_enh@result)

# A549 enhancer promoter ranking
BT549_count <- read.delim("MYC_A493_enhancer.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","MYC_A493.Tag.Count.in.1000.bp..16697286.0.Total..normalization.factor...0.60..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","MYC_A493.Tag.Count.in.1000.bp..16697286.0.Total..normalization.factor...0.60..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

A493_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
A493_prom <- as.data.frame(A493_prom@result)

A493_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
A493_enh <- as.data.frame(A493_enh@result)


# Subset adjusted p-value and pathway name for integration and join the dataset
H128_enh <- H128_enh[,c("ID","p.adjust")]
colnames(H128_enh) <- c("pathway","H128")

HT1080_enh <- HT1080_enh[,c("ID","p.adjust")]
colnames(HT1080_enh) <- c("pathway","HT1080")

HeLa_enh <- HeLa_enh[,c("ID","p.adjust")]
colnames(HeLa_enh) <- c("pathway","HeLa")

MM1S_enh <- MM1S_enh[,c("ID","p.adjust")]
colnames(MM1S_enh) <- c("pathway","MM1S")

A493_enh <- A493_enh[,c("ID","p.adjust")]
colnames(A493_enh) <- c("pathway","A493")

BT549_enh <- BT549_enh[,c("ID","p.adjust")]
colnames(BT549_enh) <- c("pathway","BT549")

MCF7_enh <- MCF7_enh[,c("ID","p.adjust")]
colnames(MCF7_enh) <- c("pathway","MCF7")

H2171_enh <- H2171_enh[,c("ID","p.adjust")]
colnames(H2171_enh) <- c("pathway","H2171")

U2OS_enh <- U2OS_enh[,c("ID","p.adjust")]
colnames(U2OS_enh) <- c("pathway","U2OS")


heatmap <- full_join(BT549_enh, MCF7_enh)
heatmap <- full_join(heatmap, H2171_enh)
heatmap <- full_join(heatmap, U2OS_enh )
heatmap <- full_join(heatmap, HT1080_enh)
heatmap <- full_join(heatmap, HeLa_enh)
heatmap <- full_join(heatmap, MM1S_enh)
heatmap <- full_join(heatmap, A493_enh)

rownames(heatmap) <- heatmap$pathway
heatmap <- heatmap[,-1]

# Select all terms which have a adj. p-value below significant
heatmap <- heatmap[heatmap$BT549 < 0.05 |
                     heatmap$MCF7 < 0.05 |
                     heatmap$H2171 < 0.05 |
                     heatmap$U2OS < 0.05 |
                     heatmap$HT1080 < 0.05 |
                     heatmap$HeLa < 0.05 |
                     heatmap$MM1S < 0.05 |
                     heatmap$A493 < 0.05 ,]
heatmap <- heatmap[complete.cases(heatmap),]

# If p-value is above 0.05 add 0.05 so it gets white in the heatmap
heatmap[heatmap > 0.05] <- 0.05 

# plot heatmaps
pheatmap(heatmap, scale = "none", color = colorRampPalette(c("#253884","white"))(200))


# Subset adjusted p-value and pathway name for integration and join the dataset
H128_prom <- H128_prom[,c("ID","p.adjust")]
colnames(H128_prom) <- c("pathway","H128")

HT1080_prom <- HT1080_prom[,c("ID","p.adjust")]
colnames(HT1080_prom) <- c("pathway","HT1080")

HeLa_prom <- HeLa_prom[,c("ID","p.adjust")]
colnames(HeLa_prom) <- c("pathway","HeLa")

MM1S_prom <- MM1S_prom[,c("ID","p.adjust")]
colnames(MM1S_prom) <- c("pathway","MM1S")

A493_prom <- A493_prom[,c("ID","p.adjust")]
colnames(A493_prom) <- c("pathway","A493")

BT549_prom <- BT549_prom[,c("ID","p.adjust")]
colnames(BT549_prom) <- c("pathway","BT549")

MCF7_prom <- MCF7_prom[,c("ID","p.adjust")]
colnames(MCF7_prom) <- c("pathway","MCF7")

H2171_prom <- H2171_prom[,c("ID","p.adjust")]
colnames(H2171_prom) <- c("pathway","H2171")

U2OS_prom <- U2OS_prom[,c("ID","p.adjust")]
colnames(U2OS_prom) <- c("pathway","U2OS")


heatmap <- full_join(BT549_prom, MCF7_prom)
heatmap <- full_join(heatmap, H2171_prom)
heatmap <- full_join(heatmap, U2OS_prom )
heatmap <- full_join(heatmap, HT1080_prom)
heatmap <- full_join(heatmap, HeLa_prom)
heatmap <- full_join(heatmap, MM1S_prom)
heatmap <- full_join(heatmap, A493_prom)


rownames(heatmap) <- heatmap$pathway
heatmap <- heatmap[,-1]

# Select all terms which have a adj. p-value below significant
heatmap <- heatmap[heatmap$BT549 < 0.05 |
                     heatmap$MCF7 < 0.05 |
                     heatmap$H2171 < 0.05 |
                     heatmap$U2OS < 0.05 |
                     heatmap$HT1080 < 0.05 |
                     heatmap$HeLa < 0.05 |
                     heatmap$MM1S < 0.05 |
                     heatmap$A493 < 0.05 ,]

heatmap <- heatmap[complete.cases(heatmap),]

# If p-value is above 0.05 add 0.05 so it gets white in the heatmap
heatmap[heatmap > 0.05] <- 0.05

# Subset adjusted p-value and pathway name for integration and join the dataset
pheatmap(heatmap, scale = "none", color = colorRampPalette(c("#B11518","white"))(200))

# Add variance
# Take gene expression data from DepMAP
Depmap <- read.csv("Expression_Public_23Q2.csv")

# Take expression from cells lines
BT549 <- Depmap[Depmap$cell_line_display_name == c("BT549"),]
MCF7 <- Depmap[Depmap$cell_line_display_name == c("MCF7"),]
U2OS <- Depmap[Depmap$cell_line_display_name == c("U2OS"),]
A549 <- Depmap[Depmap$cell_line_display_name == c("A549"),]
MM1S <- Depmap[Depmap$cell_line_display_name == c("MM1S"),]
HELA <- Depmap[Depmap$cell_line_display_name == c("HELA"),]
HT1080 <- Depmap[Depmap$cell_line_display_name == c("HT1080"),]
H2171 <- Depmap[Depmap$cell_line_display_name == c("NCIH2171"),]

Cell_line <- as.data.frame(t(rbind(BT549, MCF7,U2OS, A549, MM1S, HELA,HT1080, H2171)))
Cell_line$symbol <- colnames(Depmap)
colnames(Cell_line) <- c("BT549","MCF7","U2OS","A549","MM1S","HeLa","HT1080","H2171", "symbol")
Cell_line <- Cell_line[8:19152,]
Cell_line[, 1:8] <- sapply(Cell_line[, 1:8], as.numeric)

# Make list of promoter genes
HT1080_prom <- HT1080_prom[,c("ID","p.adjust","geneID")]
colnames(HT1080_prom) <- c("pathway","HT1080","HT1080_geneID")

HeLa_prom <- HeLa_prom[,c("ID","p.adjust","geneID")]
colnames(HeLa_prom) <- c("pathway","HeLa","HeLa_geneID")

MM1S_prom <- MM1S_prom[,c("ID","p.adjust","geneID")]
colnames(MM1S_prom) <- c("pathway","MM1S","MM1S_geneID")

A493_prom <- A493_prom[,c("ID","p.adjust","geneID")]
colnames(A493_prom) <- c("pathway","A493","A493_geneID")

BT549_prom <- BT549_prom[,c("ID","p.adjust","geneID")]
colnames(BT549_prom) <- c("pathway","BT549","BT549_geneID")

MCF7_prom <- MCF7_prom[,c("ID","p.adjust","geneID")]
colnames(MCF7_prom) <- c("pathway","MCF7","MCF7_geneID")

H2171_prom <- H2171_prom[,c("ID","p.adjust","geneID")]
colnames(H2171_prom) <- c("pathway","H2171","H2171_geneID")

U2OS_prom <- U2OS_prom[,c("ID","p.adjust","geneID")]
colnames(U2OS_prom) <- c("pathway","U2OS","U2OS_geneID")


heatmap <- full_join(BT549_prom, MCF7_prom)
heatmap <- full_join(heatmap, H2171_prom)
heatmap <- full_join(heatmap, U2OS_prom )
heatmap <- full_join(heatmap, HT1080_prom)
heatmap <- full_join(heatmap, HeLa_prom)
heatmap <- full_join(heatmap, MM1S_prom)
heatmap <- full_join(heatmap, A493_prom)


rownames(heatmap) <- heatmap$pathway
heatmap <- heatmap[,-1]


heatmap <- heatmap[heatmap$BT549 < 0.05 |
                     heatmap$MCF7 < 0.05 |
                     heatmap$H2171 < 0.05 |
                     heatmap$U2OS < 0.05 |
                     heatmap$HT1080 < 0.05 |
                     heatmap$HeLa < 0.05 |
                     heatmap$MM1S < 0.05 |
                     heatmap$A493 < 0.05 ,]

heatmap <- heatmap[1:10,c("HT1080_geneID","HeLa_geneID","MM1S_geneID","A493_geneID","BT549_geneID","MCF7_geneID","H2171_geneID","U2OS_geneID")]

# Make gene lists and subset expression
HALLMARK_MYC_TARGETS_V1 <- str_split_fixed(heatmap[c("HALLMARK_MYC_TARGETS_V1"),], ' ', 2)
HALLMARK_MYC_TARGETS_V1 <- strsplit(HALLMARK_MYC_TARGETS_V1, "/")
HALLMARK_MYC_TARGETS_V1 <- unlist(HALLMARK_MYC_TARGETS_V1)
HALLMARK_MYC_TARGETS_V1 <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_MYC_TARGETS_V1, 'SYMBOL',"ENTREZID"))

HALLMARK_MYC_TARGETS_V1 <- Cell_line[Cell_line$symbol %in% HALLMARK_MYC_TARGETS_V1$`mapIds(org.Hs.eg.db, HALLMARK_MYC_TARGETS_V1, "SYMBOL", "ENTREZID")`,]
HALLMARK_MYC_TARGETS_V1 <- as.data.frame(colMeans(HALLMARK_MYC_TARGETS_V1[,1:8]))

HALLMARK_MTORC1_SIGNALING <- str_split_fixed(heatmap[c("HALLMARK_MTORC1_SIGNALING"),], ' ', 2)
HALLMARK_MTORC1_SIGNALING <- strsplit(HALLMARK_MTORC1_SIGNALING, "/")
HALLMARK_MTORC1_SIGNALING <- unlist(HALLMARK_MTORC1_SIGNALING)
HALLMARK_MTORC1_SIGNALING <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_MTORC1_SIGNALING, 'SYMBOL',"ENTREZID"))

HALLMARK_MTORC1_SIGNALING <- Cell_line[Cell_line$symbol %in% HALLMARK_MTORC1_SIGNALING$`mapIds(org.Hs.eg.db, HALLMARK_MTORC1_SIGNALING, "SYMBOL", "ENTREZID")`,]
HALLMARK_MTORC1_SIGNALING <- as.data.frame(colMeans(HALLMARK_MTORC1_SIGNALING[,1:8]))

HALLMARK_DNA_REPAIR <- str_split_fixed(heatmap[c("HALLMARK_DNA_REPAIR"),], ' ', 2)
HALLMARK_DNA_REPAIR <- strsplit(HALLMARK_DNA_REPAIR, "/")
HALLMARK_DNA_REPAIR <- unlist(HALLMARK_DNA_REPAIR)
HALLMARK_DNA_REPAIR <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_DNA_REPAIR, 'SYMBOL',"ENTREZID"))

HALLMARK_DNA_REPAIR <- Cell_line[Cell_line$symbol %in% HALLMARK_DNA_REPAIR$`mapIds(org.Hs.eg.db, HALLMARK_DNA_REPAIR, "SYMBOL", "ENTREZID")`,]
HALLMARK_DNA_REPAIR <- as.data.frame(colMeans(HALLMARK_DNA_REPAIR[,1:8]))

HALLMARK_E2F_TARGETS <- str_split_fixed(heatmap[c("HALLMARK_E2F_TARGETS"),], ' ', 2)
HALLMARK_E2F_TARGETS <- strsplit(HALLMARK_E2F_TARGETS, "/")
HALLMARK_E2F_TARGETS <- unlist(HALLMARK_E2F_TARGETS)
HALLMARK_E2F_TARGETS <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_E2F_TARGETS, 'SYMBOL',"ENTREZID"))

HALLMARK_E2F_TARGETS <- Cell_line[Cell_line$symbol %in% HALLMARK_E2F_TARGETS$`mapIds(org.Hs.eg.db, HALLMARK_E2F_TARGETS, "SYMBOL", "ENTREZID")`,]
HALLMARK_E2F_TARGETS <- as.data.frame(colMeans(HALLMARK_E2F_TARGETS[,1:8]))

HALLMARK_G2M_CHECKPOINT <- str_split_fixed(heatmap[c("HALLMARK_G2M_CHECKPOINT"),], ' ', 2)
HALLMARK_G2M_CHECKPOINT <- strsplit(HALLMARK_G2M_CHECKPOINT, "/")
HALLMARK_G2M_CHECKPOINT <- unlist(HALLMARK_G2M_CHECKPOINT)
HALLMARK_G2M_CHECKPOINT <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_G2M_CHECKPOINT, 'SYMBOL',"ENTREZID"))

HALLMARK_G2M_CHECKPOINT <- Cell_line[Cell_line$symbol %in% HALLMARK_G2M_CHECKPOINT$`mapIds(org.Hs.eg.db, HALLMARK_G2M_CHECKPOINT, "SYMBOL", "ENTREZID")`,]
HALLMARK_G2M_CHECKPOINT <- as.data.frame(colMeans(HALLMARK_G2M_CHECKPOINT[,1:8]))

HALLMARK_MYC_TARGETS_V2 <- str_split_fixed(heatmap[c("HALLMARK_MYC_TARGETS_V2"),], ' ', 2)
HALLMARK_MYC_TARGETS_V2 <- strsplit(HALLMARK_MYC_TARGETS_V2, "/")
HALLMARK_MYC_TARGETS_V2 <- unlist(HALLMARK_MYC_TARGETS_V2)
HALLMARK_MYC_TARGETS_V2 <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_MYC_TARGETS_V2, 'SYMBOL',"ENTREZID"))

HALLMARK_MYC_TARGETS_V2 <- Cell_line[Cell_line$symbol %in% HALLMARK_MYC_TARGETS_V2$`mapIds(org.Hs.eg.db, HALLMARK_MYC_TARGETS_V2, "SYMBOL", "ENTREZID")`,]
HALLMARK_MYC_TARGETS_V2 <- as.data.frame(colMeans(HALLMARK_MYC_TARGETS_V2[,1:8]))

HALLMARK_UNFOLDED_PROTEIN_RESPONSE <- str_split_fixed(heatmap[c("HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),], ' ', 2)
HALLMARK_UNFOLDED_PROTEIN_RESPONSE <- strsplit(HALLMARK_UNFOLDED_PROTEIN_RESPONSE, "/")
HALLMARK_UNFOLDED_PROTEIN_RESPONSE <- unlist(HALLMARK_UNFOLDED_PROTEIN_RESPONSE)
HALLMARK_UNFOLDED_PROTEIN_RESPONSE <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_UNFOLDED_PROTEIN_RESPONSE, 'SYMBOL',"ENTREZID"))

HALLMARK_UNFOLDED_PROTEIN_RESPONSE <- Cell_line[Cell_line$symbol %in% HALLMARK_UNFOLDED_PROTEIN_RESPONSE$`mapIds(org.Hs.eg.db, HALLMARK_UNFOLDED_PROTEIN_RESPONSE, "SYMBOL", "ENTREZID")`,]
HALLMARK_UNFOLDED_PROTEIN_RESPONSE <- as.data.frame(colMeans(HALLMARK_UNFOLDED_PROTEIN_RESPONSE[,1:8]))

HALLMARK_OXIDATIVE_PHOSPHORYLATION <- str_split_fixed(heatmap[c("HALLMARK_OXIDATIVE_PHOSPHORYLATION"),], ' ', 2)
HALLMARK_OXIDATIVE_PHOSPHORYLATION <- strsplit(HALLMARK_OXIDATIVE_PHOSPHORYLATION, "/")
HALLMARK_OXIDATIVE_PHOSPHORYLATION <- unlist(HALLMARK_OXIDATIVE_PHOSPHORYLATION)
HALLMARK_OXIDATIVE_PHOSPHORYLATION <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_OXIDATIVE_PHOSPHORYLATION, 'SYMBOL',"ENTREZID"))

HALLMARK_OXIDATIVE_PHOSPHORYLATION <- Cell_line[Cell_line$symbol %in% HALLMARK_OXIDATIVE_PHOSPHORYLATION$`mapIds(org.Hs.eg.db, HALLMARK_OXIDATIVE_PHOSPHORYLATION, "SYMBOL", "ENTREZID")`,]
HALLMARK_OXIDATIVE_PHOSPHORYLATION <- as.data.frame(colMeans(HALLMARK_OXIDATIVE_PHOSPHORYLATION[,1:8]))

HALLMARK_PROTEIN_SECRETION <- str_split_fixed(heatmap[c("HALLMARK_PROTEIN_SECRETION"),], ' ', 2)
HALLMARK_PROTEIN_SECRETION <- strsplit(HALLMARK_PROTEIN_SECRETION, "/")
HALLMARK_PROTEIN_SECRETION <- unlist(HALLMARK_PROTEIN_SECRETION)
HALLMARK_PROTEIN_SECRETION <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_PROTEIN_SECRETION, 'SYMBOL',"ENTREZID"))

HALLMARK_PROTEIN_SECRETION <- Cell_line[Cell_line$symbol %in% HALLMARK_PROTEIN_SECRETION$`mapIds(org.Hs.eg.db, HALLMARK_PROTEIN_SECRETION, "SYMBOL", "ENTREZID")`,]
HALLMARK_PROTEIN_SECRETION <- as.data.frame(colMeans(HALLMARK_PROTEIN_SECRETION[,1:8]))

HALLMARK_MITOTIC_SPINDLE <- str_split_fixed(heatmap[c("HALLMARK_MITOTIC_SPINDLE"),], ' ', 2)
HALLMARK_MITOTIC_SPINDLE <- strsplit(HALLMARK_MITOTIC_SPINDLE, "/")
HALLMARK_MITOTIC_SPINDLE <- unlist(HALLMARK_MITOTIC_SPINDLE)
HALLMARK_MITOTIC_SPINDLE <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_MITOTIC_SPINDLE, 'SYMBOL',"ENTREZID"))

HALLMARK_MITOTIC_SPINDLE <- Cell_line[Cell_line$symbol %in% HALLMARK_MITOTIC_SPINDLE$`mapIds(org.Hs.eg.db, HALLMARK_MITOTIC_SPINDLE, "SYMBOL", "ENTREZID")`,]
HALLMARK_MITOTIC_SPINDLE <- as.data.frame(colMeans(HALLMARK_MITOTIC_SPINDLE[,1:8]))

# Merge all gene sets expression
promoter_genes <- t(data.frame(HALLMARK_MYC_TARGETS_V1,
                               HALLMARK_MTORC1_SIGNALING,
                               HALLMARK_DNA_REPAIR,
                               HALLMARK_E2F_TARGETS,
                               HALLMARK_G2M_CHECKPOINT,
                               HALLMARK_MYC_TARGETS_V2,
                               HALLMARK_UNFOLDED_PROTEIN_RESPONSE,
                               HALLMARK_OXIDATIVE_PHOSPHORYLATION,
                               HALLMARK_PROTEIN_SECRETION,
                               HALLMARK_MITOTIC_SPINDLE))

# Put back propper pathway names
rownames(promoter_genes) <- c("HALLMARK_MYC_TARGETS_V1",
                              "HALLMARK_E2F_TARGETS",
                              "HALLMARK_MTORC1_SIGNALING",
                              "HALLMARK_MYC_TARGETS_V2",
                              "HALLMARK_DNA_REPAIR",
                              "HALLMARK_G2M_CHECKPOINT",
                              "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                              "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                              "HALLMARK_PROTEIN_SECRETION",
                              "HALLMARK_MITOTIC_SPINDLE")

# Change order of plot so it fits previous order
promoter_genes <- promoter_genes[c("HALLMARK_MYC_TARGETS_V1",
                                   "HALLMARK_E2F_TARGETS",
                                   "HALLMARK_MTORC1_SIGNALING",
                                   "HALLMARK_MYC_TARGETS_V2",
                                   "HALLMARK_DNA_REPAIR",
                                   "HALLMARK_G2M_CHECKPOINT",
                                   "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                                   "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                                   "HALLMARK_PROTEIN_SECRETION",
                                   "HALLMARK_MITOTIC_SPINDLE"),]

# Valulate variance between cell lines for each term
promoter_genes <- t(rowVars(   promoter_genes)) 

# Give correct row names again
colnames(promoter_genes) <- c("HALLMARK_MYC_TARGETS_V1",
                              "HALLMARK_E2F_TARGETS",
                              "HALLMARK_MTORC1_SIGNALING",
                              "HALLMARK_MYC_TARGETS_V2",
                              "HALLMARK_DNA_REPAIR",
                              "HALLMARK_G2M_CHECKPOINT",
                              "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                              "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                              "HALLMARK_PROTEIN_SECRETION",
                              "HALLMARK_MITOTIC_SPINDLE")

# Prepare color for plotting
col_fun = colorRamp2(c(0.6, 0.3, 0), c("#3F7A13", "#3F7A13","white"))
col_fun(seq(-3, 3))

# plot heatmap
Heatmap(t(promoter_genes), name = "mat", cluster_columns = FALSE, cluster_rows = FALSE, col = col_fun)

# MSame for enhancer vairance gene expression as above
HT1080_enh <- HT1080_enh[,c("ID","p.adjust","geneID")]
colnames(HT1080_enh) <- c("pathway","HT1080","HT1080_geneID")

HeLa_enh <- HeLa_enh[,c("ID","p.adjust","geneID")]
colnames(HeLa_enh) <- c("pathway","HeLa","HeLa_geneID")

MM1S_enh <- MM1S_enh[,c("ID","p.adjust","geneID")]
colnames(MM1S_enh) <- c("pathway","MM1S","MM1S_geneID")

A493_enh <- A493_enh[,c("ID","p.adjust","geneID")]
colnames(A493_enh) <- c("pathway","A493","A493_geneID")

BT549_enh <- BT549_enh[,c("ID","p.adjust","geneID")]
colnames(BT549_enh) <- c("pathway","BT549","BT549_geneID")

MCF7_enh <- MCF7_enh[,c("ID","p.adjust","geneID")]
colnames(MCF7_enh) <- c("pathway","MCF7","MCF7_geneID")

H2171_enh <- H2171_enh[,c("ID","p.adjust","geneID")]
colnames(H2171_enh) <- c("pathway","H2171","H2171_geneID")

U2OS_enh <- U2OS_enh[,c("ID","p.adjust","geneID")]
colnames(U2OS_enh) <- c("pathway","U2OS","U2OS_geneID")


heatmap <- full_join(BT549_enh, MCF7_enh)
heatmap <- full_join(heatmap, H2171_enh)
heatmap <- full_join(heatmap, U2OS_enh )
heatmap <- full_join(heatmap, HT1080_enh)
heatmap <- full_join(heatmap, HeLa_enh)
heatmap <- full_join(heatmap, MM1S_enh)
heatmap <- full_join(heatmap, A493_enh)


rownames(heatmap) <- heatmap$pathway
heatmap <- heatmap[,-1]


heatmap <- heatmap[heatmap$BT549 < 0.05 |
                     heatmap$MCF7 < 0.05 |
                     heatmap$H2171 < 0.05 |
                     heatmap$U2OS < 0.05 |
                     heatmap$HT1080 < 0.05 |
                     heatmap$HeLa < 0.05 |
                     heatmap$MM1S < 0.05 |
                     heatmap$A493 < 0.05 ,]

heatmap <- heatmap[1:22,c("HT1080_geneID","HeLa_geneID","MM1S_geneID","A493_geneID","BT549_geneID","MCF7_geneID","H2171_geneID","U2OS_geneID")]

HALLMARK_TNFA_SIGNALING_VIA_NFKB <- str_split_fixed(heatmap[c("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),], ' ', 2)
HALLMARK_TNFA_SIGNALING_VIA_NFKB <- strsplit(HALLMARK_TNFA_SIGNALING_VIA_NFKB, "/")
HALLMARK_TNFA_SIGNALING_VIA_NFKB <- unlist(HALLMARK_TNFA_SIGNALING_VIA_NFKB)
HALLMARK_TNFA_SIGNALING_VIA_NFKB <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_TNFA_SIGNALING_VIA_NFKB, 'SYMBOL',"ENTREZID"))

HALLMARK_TNFA_SIGNALING_VIA_NFKB <- Cell_line[Cell_line$symbol %in% HALLMARK_TNFA_SIGNALING_VIA_NFKB$`mapIds(org.Hs.eg.db, HALLMARK_TNFA_SIGNALING_VIA_NFKB, "SYMBOL", "ENTREZID")`,]
HALLMARK_TNFA_SIGNALING_VIA_NFKB <- as.data.frame(colMeans(HALLMARK_TNFA_SIGNALING_VIA_NFKB[,1:8]))

HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION <- str_split_fixed(heatmap[c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),], ' ', 2)
HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION <- strsplit(HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, "/")
HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION <- unlist(HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)
HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, 'SYMBOL',"ENTREZID"))

HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION <- Cell_line[Cell_line$symbol %in% HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION$`mapIds(org.Hs.eg.db, HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, "SYMBOL", "ENTREZID")`,]
HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION <- as.data.frame(colMeans(HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION[,1:8]))

HALLMARK_HYPOXIA <- str_split_fixed(heatmap[c("HALLMARK_HYPOXIA"),], ' ', 2)
HALLMARK_HYPOXIA <- strsplit(HALLMARK_HYPOXIA, "/")
HALLMARK_HYPOXIA <- unlist(HALLMARK_HYPOXIA)
HALLMARK_HYPOXIA <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_HYPOXIA, 'SYMBOL',"ENTREZID"))

HALLMARK_HYPOXIA <- Cell_line[Cell_line$symbol %in% HALLMARK_HYPOXIA$`mapIds(org.Hs.eg.db, HALLMARK_HYPOXIA, "SYMBOL", "ENTREZID")`,]
HALLMARK_HYPOXIA <- as.data.frame(colMeans(HALLMARK_HYPOXIA[,1:8]))

HALLMARK_UV_RESPONSE_DN <- str_split_fixed(heatmap[c("HALLMARK_UV_RESPONSE_DN"),], ' ', 2)
HALLMARK_UV_RESPONSE_DN <- strsplit(HALLMARK_UV_RESPONSE_DN, "/")
HALLMARK_UV_RESPONSE_DN <- unlist(HALLMARK_UV_RESPONSE_DN)
HALLMARK_UV_RESPONSE_DN <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_UV_RESPONSE_DN, 'SYMBOL',"ENTREZID"))

HALLMARK_UV_RESPONSE_DN <- Cell_line[Cell_line$symbol %in% HALLMARK_UV_RESPONSE_DN$`mapIds(org.Hs.eg.db, HALLMARK_UV_RESPONSE_DN, "SYMBOL", "ENTREZID")`,]
HALLMARK_UV_RESPONSE_DN <- as.data.frame(colMeans(HALLMARK_UV_RESPONSE_DN[,1:8]))

HALLMARK_TGF_BETA_SIGNALING <- str_split_fixed(heatmap[c("HALLMARK_TGF_BETA_SIGNALING"),], ' ', 2)
HALLMARK_TGF_BETA_SIGNALING <- strsplit(HALLMARK_TGF_BETA_SIGNALING, "/")
HALLMARK_TGF_BETA_SIGNALING <- unlist(HALLMARK_TGF_BETA_SIGNALING)
HALLMARK_TGF_BETA_SIGNALING <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_TGF_BETA_SIGNALING, 'SYMBOL',"ENTREZID"))

HALLMARK_TGF_BETA_SIGNALING <- Cell_line[Cell_line$symbol %in% HALLMARK_TGF_BETA_SIGNALING$`mapIds(org.Hs.eg.db, HALLMARK_TGF_BETA_SIGNALING, "SYMBOL", "ENTREZID")`,]
HALLMARK_TGF_BETA_SIGNALING <- as.data.frame(colMeans(HALLMARK_TGF_BETA_SIGNALING[,1:8]))

HALLMARK_APOPTOSIS <- str_split_fixed(heatmap[c("HALLMARK_APOPTOSIS"),], ' ', 2)
HALLMARK_APOPTOSIS <- strsplit(HALLMARK_APOPTOSIS, "/")
HALLMARK_APOPTOSIS <- unlist(HALLMARK_APOPTOSIS)
HALLMARK_APOPTOSIS <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_APOPTOSIS, 'SYMBOL',"ENTREZID"))

HALLMARK_APOPTOSIS <- Cell_line[Cell_line$symbol %in% HALLMARK_APOPTOSIS$`mapIds(org.Hs.eg.db, HALLMARK_APOPTOSIS, "SYMBOL", "ENTREZID")`,]
HALLMARK_APOPTOSIS <- as.data.frame(colMeans(HALLMARK_APOPTOSIS[,1:8]))

HALLMARK_ESTROGEN_RESPONSE_LATE <- str_split_fixed(heatmap[c("HALLMARK_ESTROGEN_RESPONSE_LATE"),], ' ', 2)
HALLMARK_ESTROGEN_RESPONSE_LATE <- strsplit(HALLMARK_ESTROGEN_RESPONSE_LATE, "/")
HALLMARK_ESTROGEN_RESPONSE_LATE <- unlist(HALLMARK_ESTROGEN_RESPONSE_LATE)
HALLMARK_ESTROGEN_RESPONSE_LATE <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_ESTROGEN_RESPONSE_LATE, 'SYMBOL',"ENTREZID"))

HALLMARK_ESTROGEN_RESPONSE_LATE <- Cell_line[Cell_line$symbol %in% HALLMARK_ESTROGEN_RESPONSE_LATE$`mapIds(org.Hs.eg.db, HALLMARK_ESTROGEN_RESPONSE_LATE, "SYMBOL", "ENTREZID")`,]
HALLMARK_ESTROGEN_RESPONSE_LATE <- as.data.frame(colMeans(HALLMARK_ESTROGEN_RESPONSE_LATE[,1:8]))

HALLMARK_ANDROGEN_RESPONSE <- str_split_fixed(heatmap[c("HALLMARK_ANDROGEN_RESPONSE"),], ' ', 2)
HALLMARK_ANDROGEN_RESPONSE <- strsplit(HALLMARK_ANDROGEN_RESPONSE, "/")
HALLMARK_ANDROGEN_RESPONSE <- unlist(HALLMARK_ANDROGEN_RESPONSE)
HALLMARK_ANDROGEN_RESPONSE <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_ANDROGEN_RESPONSE, 'SYMBOL',"ENTREZID"))

HALLMARK_ANDROGEN_RESPONSE <- Cell_line[Cell_line$symbol %in% HALLMARK_ANDROGEN_RESPONSE$`mapIds(org.Hs.eg.db, HALLMARK_ANDROGEN_RESPONSE, "SYMBOL", "ENTREZID")`,]
HALLMARK_ANDROGEN_RESPONSE <- as.data.frame(colMeans(HALLMARK_ANDROGEN_RESPONSE[,1:8]))

HALLMARK_TGF_BETA_SIGNALING <- str_split_fixed(heatmap[c("HALLMARK_TGF_BETA_SIGNALING"),], ' ', 2)
HALLMARK_TGF_BETA_SIGNALING <- strsplit(HALLMARK_TGF_BETA_SIGNALING, "/")
HALLMARK_TGF_BETA_SIGNALING <- unlist(HALLMARK_TGF_BETA_SIGNALING)
HALLMARK_TGF_BETA_SIGNALING <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_TGF_BETA_SIGNALING, 'SYMBOL',"ENTREZID"))

HALLMARK_TGF_BETA_SIGNALING <- Cell_line[Cell_line$symbol %in% HALLMARK_TGF_BETA_SIGNALING$`mapIds(org.Hs.eg.db, HALLMARK_TGF_BETA_SIGNALING, "SYMBOL", "ENTREZID")`,]
HALLMARK_TGF_BETA_SIGNALING <- as.data.frame(colMeans(HALLMARK_TGF_BETA_SIGNALING[,1:8]))

HALLMARK_ESTROGEN_RESPONSE_EARLY <- str_split_fixed(heatmap[c("HALLMARK_ESTROGEN_RESPONSE_EARLY"),], ' ', 2)
HALLMARK_ESTROGEN_RESPONSE_EARLY <- strsplit(HALLMARK_ESTROGEN_RESPONSE_EARLY, "/")
HALLMARK_ESTROGEN_RESPONSE_EARLY <- unlist(HALLMARK_ESTROGEN_RESPONSE_EARLY)
HALLMARK_ESTROGEN_RESPONSE_EARLY <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_ESTROGEN_RESPONSE_EARLY, 'SYMBOL',"ENTREZID"))

HALLMARK_ESTROGEN_RESPONSE_EARLY <- Cell_line[Cell_line$symbol %in% HALLMARK_ESTROGEN_RESPONSE_EARLY$`mapIds(org.Hs.eg.db, HALLMARK_ESTROGEN_RESPONSE_EARLY, "SYMBOL", "ENTREZID")`,]
HALLMARK_ESTROGEN_RESPONSE_EARLY <- as.data.frame(colMeans(HALLMARK_ESTROGEN_RESPONSE_EARLY[,1:8]))


HALLMARK_APICAL_JUNCTION <- str_split_fixed(heatmap[c("HALLMARK_APICAL_JUNCTION"),], ' ', 2)
HALLMARK_APICAL_JUNCTION <- strsplit(HALLMARK_APICAL_JUNCTION, "/")
HALLMARK_APICAL_JUNCTION <- unlist(HALLMARK_APICAL_JUNCTION)
HALLMARK_APICAL_JUNCTION <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_APICAL_JUNCTION, 'SYMBOL',"ENTREZID"))

HALLMARK_APICAL_JUNCTION <- Cell_line[Cell_line$symbol %in% HALLMARK_APICAL_JUNCTION$`mapIds(org.Hs.eg.db, HALLMARK_APICAL_JUNCTION, "SYMBOL", "ENTREZID")`,]
HALLMARK_APICAL_JUNCTION <- as.data.frame(colMeans(HALLMARK_APICAL_JUNCTION[,1:8]))


HALLMARK_IL2_STAT5_SIGNALING <- str_split_fixed(heatmap[c("HALLMARK_IL2_STAT5_SIGNALING"),], ' ', 2)
HALLMARK_IL2_STAT5_SIGNALING <- strsplit(HALLMARK_IL2_STAT5_SIGNALING, "/")
HALLMARK_IL2_STAT5_SIGNALING <- unlist(HALLMARK_IL2_STAT5_SIGNALING)
HALLMARK_IL2_STAT5_SIGNALING <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_IL2_STAT5_SIGNALING, 'SYMBOL',"ENTREZID"))

HALLMARK_IL2_STAT5_SIGNALING <- Cell_line[Cell_line$symbol %in% HALLMARK_IL2_STAT5_SIGNALING$`mapIds(org.Hs.eg.db, HALLMARK_IL2_STAT5_SIGNALING, "SYMBOL", "ENTREZID")`,]
HALLMARK_IL2_STAT5_SIGNALING <- as.data.frame(colMeans(HALLMARK_IL2_STAT5_SIGNALING[,1:8]))


HALLMARK_MITOTIC_SPINDLE <- str_split_fixed(heatmap[c("HALLMARK_MITOTIC_SPINDLE"),], ' ', 2)
HALLMARK_MITOTIC_SPINDLE <- strsplit(HALLMARK_MITOTIC_SPINDLE, "/")
HALLMARK_MITOTIC_SPINDLE <- unlist(HALLMARK_MITOTIC_SPINDLE)
HALLMARK_MITOTIC_SPINDLE <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_MITOTIC_SPINDLE, 'SYMBOL',"ENTREZID"))

HALLMARK_MITOTIC_SPINDLE <- Cell_line[Cell_line$symbol %in% HALLMARK_MITOTIC_SPINDLE$`mapIds(org.Hs.eg.db, HALLMARK_MITOTIC_SPINDLE, "SYMBOL", "ENTREZID")`,]
HALLMARK_MITOTIC_SPINDLE <- as.data.frame(colMeans(HALLMARK_MITOTIC_SPINDLE[,1:8]))


HALLMARK_P53_PATHWAY <- str_split_fixed(heatmap[c("HALLMARK_P53_PATHWAY"),], ' ', 2)
HALLMARK_P53_PATHWAY <- strsplit(HALLMARK_P53_PATHWAY, "/")
HALLMARK_P53_PATHWAY <- unlist(HALLMARK_P53_PATHWAY)
HALLMARK_P53_PATHWAY <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_P53_PATHWAY, 'SYMBOL',"ENTREZID"))

HALLMARK_P53_PATHWAY <- Cell_line[Cell_line$symbol %in% HALLMARK_P53_PATHWAY$`mapIds(org.Hs.eg.db, HALLMARK_P53_PATHWAY, "SYMBOL", "ENTREZID")`,]
HALLMARK_P53_PATHWAY <- as.data.frame(colMeans(HALLMARK_P53_PATHWAY[,1:8]))

HALLMARK_ALLOGRAFT_REJECTION <- str_split_fixed(heatmap[c("HALLMARK_ALLOGRAFT_REJECTION"),], ' ', 2)
HALLMARK_ALLOGRAFT_REJECTION <- strsplit(HALLMARK_ALLOGRAFT_REJECTION, "/")
HALLMARK_ALLOGRAFT_REJECTION <- unlist(HALLMARK_ALLOGRAFT_REJECTION)
HALLMARK_ALLOGRAFT_REJECTION <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_ALLOGRAFT_REJECTION, 'SYMBOL',"ENTREZID"))

HALLMARK_ALLOGRAFT_REJECTION <- Cell_line[Cell_line$symbol %in% HALLMARK_ALLOGRAFT_REJECTION$`mapIds(org.Hs.eg.db, HALLMARK_ALLOGRAFT_REJECTION, "SYMBOL", "ENTREZID")`,]
HALLMARK_ALLOGRAFT_REJECTION <- as.data.frame(colMeans(HALLMARK_ALLOGRAFT_REJECTION[,1:8]))


HALLMARK_MTORC1_SIGNALING <- str_split_fixed(heatmap[c("HALLMARK_MTORC1_SIGNALING"),], ' ', 2)
HALLMARK_MTORC1_SIGNALING <- strsplit(HALLMARK_MTORC1_SIGNALING, "/")
HALLMARK_MTORC1_SIGNALING <- unlist(HALLMARK_MTORC1_SIGNALING)
HALLMARK_MTORC1_SIGNALING <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_MTORC1_SIGNALING, 'SYMBOL',"ENTREZID"))

HALLMARK_MTORC1_SIGNALING <- Cell_line[Cell_line$symbol %in% HALLMARK_MTORC1_SIGNALING$`mapIds(org.Hs.eg.db, HALLMARK_MTORC1_SIGNALING, "SYMBOL", "ENTREZID")`,]
HALLMARK_MTORC1_SIGNALING <- as.data.frame(colMeans(HALLMARK_MTORC1_SIGNALING[,1:8]))


HALLMARK_MYOGENESIS <- str_split_fixed(heatmap[c("HALLMARK_MYOGENESIS"),], ' ', 2)
HALLMARK_MYOGENESIS <- strsplit(HALLMARK_MYOGENESIS, "/")
HALLMARK_MYOGENESIS <- unlist(HALLMARK_MYOGENESIS)
HALLMARK_MYOGENESIS <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_MYOGENESIS, 'SYMBOL',"ENTREZID"))

HALLMARK_MYOGENESIS <- Cell_line[Cell_line$symbol %in% HALLMARK_MYOGENESIS$`mapIds(org.Hs.eg.db, HALLMARK_MYOGENESIS, "SYMBOL", "ENTREZID")`,]
HALLMARK_MYOGENESIS <- as.data.frame(colMeans(HALLMARK_MYOGENESIS[,1:8]))

HALLMARK_UV_RESPONSE_UP <- str_split_fixed(heatmap[c("HALLMARK_UV_RESPONSE_UP"),], ' ', 2)
HALLMARK_UV_RESPONSE_UP <- strsplit(HALLMARK_UV_RESPONSE_UP, "/")
HALLMARK_UV_RESPONSE_UP <- unlist(HALLMARK_UV_RESPONSE_UP)
HALLMARK_UV_RESPONSE_UP <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_UV_RESPONSE_UP, 'SYMBOL',"ENTREZID"))

HALLMARK_UV_RESPONSE_UP <- Cell_line[Cell_line$symbol %in% HALLMARK_UV_RESPONSE_UP$`mapIds(org.Hs.eg.db, HALLMARK_UV_RESPONSE_UP, "SYMBOL", "ENTREZID")`,]
HALLMARK_UV_RESPONSE_UP <- as.data.frame(colMeans(HALLMARK_UV_RESPONSE_UP[,1:8]))

HALLMARK_INTERFERON_GAMMA_RESPONSE <- str_split_fixed(heatmap[c("HALLMARK_INTERFERON_GAMMA_RESPONSE"),], ' ', 2)
HALLMARK_INTERFERON_GAMMA_RESPONSE <- strsplit(HALLMARK_INTERFERON_GAMMA_RESPONSE, "/")
HALLMARK_INTERFERON_GAMMA_RESPONSE <- unlist(HALLMARK_INTERFERON_GAMMA_RESPONSE)
HALLMARK_INTERFERON_GAMMA_RESPONSE <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_INTERFERON_GAMMA_RESPONSE, 'SYMBOL',"ENTREZID"))

HALLMARK_INTERFERON_GAMMA_RESPONSE <- Cell_line[Cell_line$symbol %in% HALLMARK_INTERFERON_GAMMA_RESPONSE$`mapIds(org.Hs.eg.db, HALLMARK_INTERFERON_GAMMA_RESPONSE, "SYMBOL", "ENTREZID")`,]
HALLMARK_INTERFERON_GAMMA_RESPONSE <- as.data.frame(colMeans(HALLMARK_INTERFERON_GAMMA_RESPONSE[,1:8]))

HALLMARK_UNFOLDED_PROTEIN_RESPONSE <- str_split_fixed(heatmap[c("HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),], ' ', 2)
HALLMARK_UNFOLDED_PROTEIN_RESPONSE <- strsplit(HALLMARK_UNFOLDED_PROTEIN_RESPONSE, "/")
HALLMARK_UNFOLDED_PROTEIN_RESPONSE <- unlist(HALLMARK_UNFOLDED_PROTEIN_RESPONSE)
HALLMARK_UNFOLDED_PROTEIN_RESPONSE <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_UNFOLDED_PROTEIN_RESPONSE, 'SYMBOL',"ENTREZID"))

HALLMARK_UNFOLDED_PROTEIN_RESPONSE <- Cell_line[Cell_line$symbol %in% HALLMARK_UNFOLDED_PROTEIN_RESPONSE$`mapIds(org.Hs.eg.db, HALLMARK_UNFOLDED_PROTEIN_RESPONSE, "SYMBOL", "ENTREZID")`,]
HALLMARK_UNFOLDED_PROTEIN_RESPONSE <- as.data.frame(colMeans(HALLMARK_UNFOLDED_PROTEIN_RESPONSE[,1:8]))

HALLMARK_HEDGEHOG_SIGNALING <- str_split_fixed(heatmap[c("HALLMARK_HEDGEHOG_SIGNALING"),], ' ', 2)
HALLMARK_HEDGEHOG_SIGNALING <- strsplit(HALLMARK_HEDGEHOG_SIGNALING, "/")
HALLMARK_HEDGEHOG_SIGNALING <- unlist(HALLMARK_HEDGEHOG_SIGNALING)
HALLMARK_HEDGEHOG_SIGNALING <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_HEDGEHOG_SIGNALING, 'SYMBOL',"ENTREZID"))

HALLMARK_HEDGEHOG_SIGNALING <- Cell_line[Cell_line$symbol %in% HALLMARK_HEDGEHOG_SIGNALING$`mapIds(org.Hs.eg.db, HALLMARK_HEDGEHOG_SIGNALING, "SYMBOL", "ENTREZID")`,]
HALLMARK_HEDGEHOG_SIGNALING <- as.data.frame(colMeans(HALLMARK_HEDGEHOG_SIGNALING[,1:8]))

HALLMARK_G2M_CHECKPOINT <- str_split_fixed(heatmap[c("HALLMARK_G2M_CHECKPOINT"),], ' ', 2)
HALLMARK_G2M_CHECKPOINT <- strsplit(HALLMARK_G2M_CHECKPOINT, "/")
HALLMARK_G2M_CHECKPOINT <- unlist(HALLMARK_G2M_CHECKPOINT)
HALLMARK_G2M_CHECKPOINT <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_G2M_CHECKPOINT, 'SYMBOL',"ENTREZID"))

HALLMARK_G2M_CHECKPOINT <- Cell_line[Cell_line$symbol %in% HALLMARK_G2M_CHECKPOINT$`mapIds(org.Hs.eg.db, HALLMARK_G2M_CHECKPOINT, "SYMBOL", "ENTREZID")`,]
HALLMARK_G2M_CHECKPOINT <- as.data.frame(colMeans(HALLMARK_G2M_CHECKPOINT[,1:8]))

HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY <- str_split_fixed(heatmap[c("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"),], ' ', 2)
HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY <- strsplit(HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY, "/")
HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY <- unlist(HALLMARK_TGF_BETA_SIGNALING)
HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY <- as.data.frame(mapIds(org.Hs.eg.db, HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY, 'SYMBOL',"ENTREZID"))

HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY <- Cell_line[Cell_line$symbol %in% HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY$`mapIds(org.Hs.eg.db, HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY, "SYMBOL", "ENTREZID")`,]
HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY <- as.data.frame(colMeans(HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY[,1:8]))

enhancer_genes <- t(data.frame(HALLMARK_TNFA_SIGNALING_VIA_NFKB,
                               HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
                               HALLMARK_HYPOXIA,
                               HALLMARK_UV_RESPONSE_DN,
                               HALLMARK_TGF_BETA_SIGNALING,
                               HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY,
                               HALLMARK_APOPTOSIS,
                               HALLMARK_ESTROGEN_RESPONSE_LATE,
                               HALLMARK_ANDROGEN_RESPONSE,
                               HALLMARK_ESTROGEN_RESPONSE_EARLY,
                               HALLMARK_APICAL_JUNCTION,
                               HALLMARK_IL2_STAT5_SIGNALING,
                               HALLMARK_MITOTIC_SPINDLE,
                               HALLMARK_P53_PATHWAY,
                               HALLMARK_ALLOGRAFT_REJECTION,
                               HALLMARK_MTORC1_SIGNALING,
                               HALLMARK_MYOGENESIS,
                               HALLMARK_UV_RESPONSE_UP,
                               HALLMARK_INTERFERON_GAMMA_RESPONSE,
                               HALLMARK_UNFOLDED_PROTEIN_RESPONSE,
                               HALLMARK_HEDGEHOG_SIGNALING,
                               HALLMARK_G2M_CHECKPOINT))

rownames(enhancer_genes) <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                              "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                              "HALLMARK_HYPOXIA",
                              "HALLMARK_UV_RESPONSE_DN",
                              "HALLMARK_TGF_BETA_SIGNALING",
                              "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                              "HALLMARK_APOPTOSIS",
                              "HALLMARK_ESTROGEN_RESPONSE_LATE",
                              "HALLMARK_ANDROGEN_RESPONSE",
                              "HALLMARK_ESTROGEN_RESPONSE_EARLY",
                              "HALLMARK_APICAL_JUNCTION",
                              "HALLMARK_IL2_STAT5_SIGNALING",
                              "HALLMARK_MITOTIC_SPINDLE",
                              "HALLMARK_P53_PATHWAY",
                              "HALLMARK_ALLOGRAFT_REJECTION",
                              "HALLMARK_MTORC1_SIGNALING",
                              "HALLMARK_MYOGENESIS",
                              "HALLMARK_UV_RESPONSE_UP",
                              "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                              "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                              "HALLMARK_HEDGEHOG_SIGNALING",
                              "HALLMARK_G2M_CHECKPOINT")




enhancer_genes <- enhancer_genes[c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                                   "HALLMARK_UV_RESPONSE_DN",
                                   "HALLMARK_ESTROGEN_RESPONSE_EARLY",
                                   "HALLMARK_ESTROGEN_RESPONSE_LATE",
                                   "HALLMARK_ANDROGEN_RESPONSE",
                                   "HALLMARK_TGF_BETA_SIGNALING",
                                   "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                   "HALLMARK_HYPOXIA",
                                   "HALLMARK_IL2_STAT5_SIGNALING",
                                   "HALLMARK_ALLOGRAFT_REJECTION",
                                   "HALLMARK_P53_PATHWAY",
                                   "HALLMARK_APICAL_JUNCTION",
                                   "HALLMARK_MYOGENESIS",
                                   "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                                   "HALLMARK_MITOTIC_SPINDLE",
                                   "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                                   "HALLMARK_MTORC1_SIGNALING",
                                   "HALLMARK_UV_RESPONSE_UP",
                                   "HALLMARK_APOPTOSIS",
                                   "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                                   "HALLMARK_HEDGEHOG_SIGNALING",
                                   "HALLMARK_G2M_CHECKPOINT"),]




enhancer_var <- t(rowVars(enhancer_genes)) 
colnames(enhancer_var) <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                            "HALLMARK_UV_RESPONSE_DN",
                            "HALLMARK_ESTROGEN_RESPONSE_EARLY",
                            "HALLMARK_ESTROGEN_RESPONSE_LATE",
                            "HALLMARK_ANDROGEN_RESPONSE",
                            "HALLMARK_TGF_BETA_SIGNALING",
                            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                            "HALLMARK_HYPOXIA",
                            "HALLMARK_IL2_STAT5_SIGNALING",
                            "HALLMARK_ALLOGRAFT_REJECTION",
                            "HALLMARK_P53_PATHWAY",
                            "HALLMARK_APICAL_JUNCTION",
                            "HALLMARK_MYOGENESIS",
                            "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                            "HALLMARK_MITOTIC_SPINDLE",
                            "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                            "HALLMARK_MTORC1_SIGNALING",
                            "HALLMARK_UV_RESPONSE_UP",
                            "HALLMARK_APOPTOSIS",
                            "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                            "HALLMARK_HEDGEHOG_SIGNALING",
                            "HALLMARK_G2M_CHECKPOINT")

col_fun = colorRamp2(c(0.6, 0.3, 0), c("#3F7A13", "#3F7A13","white"))
col_fun(seq(-3, 3))

Heatmap(t(enhancer_var), name = "mat", cluster_columns = FALSE, cluster_rows = FALSE, col = col_fun)




################# Figure 2F ################## 
# Run this if download does not work. The Metabrix gene expression matrix can be downloaded from cBioportal webpage. 
mb_exp2 <- read.delim("metabrix_expression.txt", sep = " ")

#Import metastasis data from Rueda et al 2019. This is clinical data for the metabric cohort downloaded as suppl table S6 from Rueda et al, Nature 2019. The death column is the same as from cbioportal as expected. Info on each column can be found in 41586_2019_1007_MOESM10_ESM_Metabric_variableDescription.xls.
relapse<-read.delim("NIHMS1520488-supplement-Supp_Table_6.txt",h=T) 

# Divide genes into quantiles and subset enhancer genes which has been connected to a genes (this plot4 is from figure 2A where we have ranked genes). Only genes which has are activated by MYC is included in this analysis
plot4_1 <- plot4[plot4$Log2FC_MYC_vs_Crtl < -0.5,]
quantile(plot4_1$count, probs = c(0.3,0.7))

quant_1 <- plot4_1[plot4_1$count > 86.100,]
quant_4 <- plot4_1[plot4_1$count < -33.782 & plot4_1$ABC.Score > 0,]


##### BT549 specific Enhancer regulated genes
# Subset genes from the metabric dataset and add clinical information 
gene_info<-as.data.frame(t(mb_exp2[rownames(mb_exp2) %in% quant_4$symbol,])) #Get gene expression for selected genes
gene_info$METABRIC.ID<-rownames(gene_info)
relapse2<-merge(relapse,gene_info,by="METABRIC.ID")

# Calculate the mean of genes in each patients
relapse2$Mean <- colMeans((t(relapse2[,36:ncol(relapse2)]))) 

# Subset patients with IDC negative for ER, PR and HER2
km<-relapse2 %>% filter(Histological.Type=="IDC", ER.Expr=="-",PR.Expr=="-",Her2.Expr=="-")

# Divide patients into high and low expression of the target enhancer genes
km$type<-"expression"
km[km$Mean>quantile(km$Mean, probs=0.5),"type"]<-"High" 
km[km$Mean<quantile(km$Mean, probs=0.5),"type"]<-"Low" 

km<-km[km$type!="expression" ,]

km_fit <- survfit(Surv(TDR, DR) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)

# Subset patients with IDC positive for ER and negative for HER HER2
km<-relapse2 %>% filter(Histological.Type=="IDC", ER.Expr=="+", Her2.Expr=="-")

# Divide patients into high and low expression of the target enhancer genes
km$type<-"expression"
km[km$Mean>quantile(km$Mean, probs=0.5),"type"]<-"High" 
km[km$Mean<quantile(km$Mean, probs=0.5),"type"]<-"Low" 

km<-km[km$type!="expression" ,]

km_fit <- survfit(Surv(TDR, DR) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)


##### BT549 promoters regulated genes
# Subset genes from the metabric dataset and add clinical information 
gene_info<-as.data.frame(t(mb_exp2[rownames(mb_exp2) %in% quant_1$symbol,])) #Get gene expression for selected genes
gene_info$METABRIC.ID<-rownames(gene_info)
relapse2<-merge(relapse,gene_info,by="METABRIC.ID")

# Calculate the mean of genes in each patients
relapse2$Mean <- colMeans((t(relapse2[,36:ncol(relapse2)]))) 

# Subset patients with IDC negative for ER, PR and HER2
km<-relapse2 %>% filter(Histological.Type=="IDC", ER.Expr=="-",PR.Expr=="-",Her2.Expr=="-")

# Divide patients into high and low expression of the target enhancer genes
km$type<-"expression"
km[km$Mean>quantile(km$Mean, probs=0.5),"type"]<-"High" 
km[km$Mean<quantile(km$Mean, probs=0.5),"type"]<-"Low" 

km<-km[km$type!="expression" ,]

km_fit <- survfit(Surv(TDR, DR) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)

# Subset patients with IDC positive for ER and negative for HER HER2
km<-relapse2 %>% filter(Histological.Type=="IDC", ER.Expr=="+", Her2.Expr=="-")
km <- km[complete.cases(km$Mean),]

# Divide patients into high and low expression of the target enhancer genes
km$type<-"expression"
km[km$Mean>quantile(km$Mean, probs=0.5),"type"]<-"High" 
km[km$Mean<quantile(km$Mean, probs=0.5),"type"]<-"Low" 

km<-km[km$type!="expression" ,]

km_fit <- survfit(Surv(TDR, DR) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)






################# Figure 2G ################## 
# Make the ranking for MYC ChIP-seq in MCF7 as for BT549 in figure 2A
MCF7_count <- read.delim("MYC_MCF7_Consensus.count.txt")

MCF7_count$Ann <- sub(" .*","",MCF7_count$Annotation)
MCF7_anno_Enhancer <- MCF7_count[MCF7_count$Ann == "Intergenic" | MCF7_count$Ann == "intron",]
MCF7_anno_Enhancer <- MCF7_anno_Enhancer[(MCF7_anno_Enhancer$Distance.to.TSS < 50000 & MCF7_anno_Enhancer$Detailed.Annotation > 3000) & MCF7_anno_Enhancer$Distance.to.TSS > -50000,]
MCF7_anno_Promoter <- MCF7_count[MCF7_count$Ann == "promoter-TSS",]
MCF7_anno_Promoter <- MCF7_anno_Promoter[!duplicated(MCF7_anno_Promoter$Gene.Name),]


MCF7_anno_Enhancer_strengh <- MCF7_anno_Enhancer[,c("Gene.Name","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_MC7..Tag.Count.in.1000.bp..42607752.5.Total..normalization.factor...0.23..effective.total...10000000.")]
MCF7_anno_Enhancer_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> MCF7_anno_Enhancer_strengh
colnames(MCF7_anno_Enhancer_strengh) <- c("symbol","Count_enh")

MCF7_anno_promoter_strengh <- MCF7_anno_Promoter[,c("Gene.Name","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_MC7..Tag.Count.in.1000.bp..42607752.5.Total..normalization.factor...0.23..effective.total...10000000.")]
MCF7_anno_promoter_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> MCF7_anno_promoter_strengh
colnames(MCF7_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(MCF7_anno_Enhancer_strengh, MCF7_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh

# Add gene expression information to the genes with nearby MYC binding if the gene is significant regulated upon MYCMI6 treatment in MCF7. The dataset is downloaded from GSE197062ID
plot2 <- read.table("DE_analysis_MYCMI6_MCF7.txt")
plot2 <- plot2[plot2$padj.c1 < .05,]
plot2 <- plot2[,c("symbol", "log2FC.c1")]
plot3 <- left_join(plot, plot2)
plot3[is.na(plot3)] = 0
plot3 <- plot3[!duplicated(plot3$symbol),]

# Load in ABC analysis which has been performed with HI-C from GSE66733
ABC_MCF7 <- read.delim("MCF7_MED_EnhancerPredictionsFull.txt")
ABC_MCF7_MYC <- read.delim("EnhancerPredictionsFull_MYC_binding_MCF7.txt")

ABC_merged <- merge(ABC_MCF7, ABC_MCF7_MYC, by.x = "name", by.y = "PeakID..cmd.annotatePeaks.pl.MCF7_Predictions_MED1.MCF7_MED_EnhancerPredictionsFull.txt.hg38..d....CellLine_MYC_MED1_ChIP.tagdir.Combined_MYC_MC7...cpu.10.")

ABC_merged <- ABC_merged[ABC_merged$...CellLine_MYC_MED1_ChIP.tagdir.Combined_MYC_MC7..Tag.Count.in.given.bp..42607752.5.Total..normalization.factor...0.23..effective.total...10000000. > 30,]

ABC_plot <- ABC_merged[,c("TargetGene","ABC.Score","powerlaw_contact")]
ABC_plot %>% group_by(TargetGene) %>% summarise_all(funs(sum)) -> ABC_plot_1
colnames(ABC_plot_1) <- c("symbol","ABC.Score","HIC")


plot4 <- left_join(plot3, ABC_plot_1)
plot4[is.na(plot4)] = 0

# Outcome for MCF7 identified genes
# Divide genes into quantiles and subset enhancer genes which has been connected to a genes (this plot4 is from figure 2A where we have ranked genes). Only genes which has are activated by MYC is included in this analysis
plot4_1 <- plot4[plot4$log2FC.c1 < -0.5,]
quantile(plot4_1$count, probs = c(0.3,0.7))

quant_1 <- plot4_1[plot4_1$count > 35.67,]
quant_4 <- plot4_1[plot4_1$count < -69.83 & plot4_1$ABC.Score > 0,]


##### MCF7 specific Enhancer regulated genes
# Subset genes from the metabric dataset and add clinical information 
gene_info<-as.data.frame(t(mb_exp2[rownames(mb_exp2) %in% quant_4$symbol,])) #Get gene expression for selected genes
gene_info$METABRIC.ID<-rownames(gene_info)
relapse2<-merge(relapse,gene_info,by="METABRIC.ID")

# Calculate the mean of genes in each patients
relapse2$Mean <- colMeans((t(relapse2[,36:ncol(relapse2)]))) 

# Subset patients with IDC negative for ER, PR and HER2
km<-relapse2 %>% filter(Histological.Type=="IDC", ER.Expr=="-",PR.Expr=="-",Her2.Expr=="-")

# Divide patients into high and low expression of the target enhancer genes
km$type<-"expression"
km[km$Mean>quantile(km$Mean, probs=0.5),"type"]<-"High" 
km[km$Mean<quantile(km$Mean, probs=0.5),"type"]<-"Low" 

km<-km[km$type!="expression" ,]

km_fit <- survfit(Surv(TDR, DR) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)

# Subset patients with IDC positive for ER and negative for HER HER2
km<-relapse2 %>% filter(Histological.Type=="IDC", ER.Expr=="+", Her2.Expr=="-")

# Divide patients into high and low expression of the target enhancer genes
km$type<-"expression"
km[km$Mean>quantile(km$Mean, probs=0.5),"type"]<-"High" 
km[km$Mean<quantile(km$Mean, probs=0.5),"type"]<-"Low" 

km<-km[km$type!="expression" ,]

km_fit <- survfit(Surv(TDR, DR) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)



##### MCF7 promoters regulated genes
# Subset genes from the metabric dataset and add clinical information 
gene_info<-as.data.frame(t(mb_exp2[rownames(mb_exp2) %in% quant_1$symbol,])) #Get gene expression for selected genes
gene_info$METABRIC.ID<-rownames(gene_info)
relapse2<-merge(relapse,gene_info,by="METABRIC.ID")

# Calculate the mean of genes in each patients
relapse2$Mean <- colMeans((t(relapse2[,36:ncol(relapse2)]))) 

# Subset patients with IDC negative for ER, PR and HER2
km<-relapse2 %>% filter(Histological.Type=="IDC", ER.Expr=="-",PR.Expr=="-",Her2.Expr=="-")

# Divide patients into high and low expression of the target enhancer genes
km$type<-"expression"
km[km$Mean>quantile(km$Mean, probs=0.5),"type"]<-"High" 
km[km$Mean<quantile(km$Mean, probs=0.5),"type"]<-"Low" 

km<-km[km$type!="expression" ,]

km_fit <- survfit(Surv(TDR, DR) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)

# Subset patients with IDC positive for ER and negative for HER HER2
km<-relapse2 %>% filter(Histological.Type=="IDC", ER.Expr=="+", Her2.Expr=="-")
km <- km[complete.cases(km$Mean),]

# Divide patients into high and low expression of the target enhancer genes
km$type<-"expression"
km[km$Mean>quantile(km$Mean, probs=0.5),"type"]<-"High" 
km[km$Mean<quantile(km$Mean, probs=0.5),"type"]<-"Low" 

km<-km[km$type!="expression" ,]

km_fit <- survfit(Surv(TDR, DR) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)










































################# Supp 2A ################## 
BT549_count <- read.delim("MYC_BT549_Consensus.count.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Gene.Name","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_BT549..Tag.Count.in.1000.bp..40939315.0.Total..normalization.factor...0.24..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Gene.Name","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_BT549..Tag.Count.in.1000.bp..40939315.0.Total..normalization.factor...0.24..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh

plot2 <- read.table("DE_analysis_KO_BT549.txt")
plot2 <- plot2[plot2$padj_MYC_vs_Crtl < .05,]
plot2 <- plot2[,c("symbol", "Log2FC_MYC_vs_Crtl")]
plot3 <- left_join(plot, plot2)
plot3[is.na(plot3)] = 0
plot3 <- plot3[!duplicated(plot3$symbol),]

plot4 <- plot3[1:16350,2:5]
plot4 <- plot4[order(plot4$count),]
plot5 <- rowsum(plot4,rep(1:327,each=50))

### Add enhancer connections
gene_enhancer <- read.delim("MYCall_enhancer_gene.txt", header = F)
gene_enhancer <- gene_enhancer[,c("V9","V1")]
gene_enhancer$V1 <- rep(1)
colnames(gene_enhancer) <- c("symbol","test")
test <- left_join(plot3, gene_enhancer)
test[is.na(test)] <- 0
test$scale <- rep("Connection")
test <- test[!duplicated(test$symbol),]
test <- test[test$Log2FC_MYC_vs_Crtl < 0 | test$Log2FC_MYC_vs_Crtl > 0, ]

ggplot(data= test,mapping=aes(x=rank(-count),y=test)) + 
  geom_smooth() + theme(axis.title = element_text(size = 20),
                        axis.text = element_text(size = 20, colour = "black"),
                        panel.background = element_blank(),
                        panel.grid = element_blank(),
                        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + xlab("Ranked Promoter - Enhancer dominace") + ylab("Density of Enh-Pro connections") + coord_cartesian(ylim=c(0,0.3))

################# Supp 2C ################## 
# Get gene set data 
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

# Subset genes which are significant up or downregulated with MYC siRNA knockdown in BT549 cells
plot_sig <- plot3[abs(plot3$Log2FC_MYC_vs_Crtl) > 0,]

plot_prom <- plot_sig[plot_sig$count>quantile(plot_sig$count, probs=0.75),"symbol"]
plot_enh <- plot_sig[-plot_sig$count>quantile(-plot_sig$count, probs=0.75),"symbol"]

plot_prom <- as.data.frame(mapIds(org.Hs.eg.db, plot_prom$symbol, "ENTREZID",'SYMBOL'))
plot_enh <- as.data.frame(mapIds(org.Hs.eg.db, plot_enh$symbol, "ENTREZID",'SYMBOL'))

BT549_prom <- enricher(plot_prom$`mapIds(org.Hs.eg.db, plot_prom$symbol, "ENTREZID", "SYMBOL")`, TERM2GENE=m_t2g)
BT549_prom <- as.data.frame(BT549_prom@result)
BT549_prom_sig <- BT549_prom[BT549_prom$p.adjust < 0.05,]

ggplot(data=BT549_prom_sig, aes(x=reorder(ID, -log10(p.adjust)), y=-log10(p.adjust))) +geom_bar(stat="identity", position = "dodge") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + coord_flip() + xlab(" ") + ylab("-log10(p.adjust)")


BT549_enh <- enricher(plot_enh$`mapIds(org.Hs.eg.db, plot_enh$symbol, "ENTREZID", "SYMBOL")`, TERM2GENE=m_t2g)
BT549_enh <- as.data.frame(BT549_enh@result)
BT549_enh_sig <- BT549_enh[BT549_enh$p.adjust < 0.05,]

ggplot(data=BT549_enh_sig, aes(x=reorder(ID, -log10(p.adjust)), y=-log10(p.adjust))) +geom_bar(stat="identity", position = "dodge") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + scale_fill_viridis_c(option="E") + coord_flip() + xlab(" ") + ylab("-log10(p.adjust)")


################# Supp 2D ##################
# Make the ranking for MYC ChIP-seq in MCF7 as for BT549 in figure 2A
MCF7_count <- read.delim("MYC_MCF7_Consensus.count.txt")

MCF7_count$Ann <- sub(" .*","",MCF7_count$Annotation)
MCF7_anno_Enhancer <- MCF7_count[MCF7_count$Ann == "Intergenic" | MCF7_count$Ann == "intron",]
MCF7_anno_Enhancer <- MCF7_anno_Enhancer[(MCF7_anno_Enhancer$Distance.to.TSS < 50000 & MCF7_anno_Enhancer$Detailed.Annotation > 3000) & MCF7_anno_Enhancer$Distance.to.TSS > -50000,]
MCF7_anno_Promoter <- MCF7_count[MCF7_count$Ann == "promoter-TSS",]
MCF7_anno_Promoter <- MCF7_anno_Promoter[!duplicated(MCF7_anno_Promoter$Gene.Name),]


MCF7_anno_Enhancer_strengh <- MCF7_anno_Enhancer[,c("Gene.Name","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_MC7..Tag.Count.in.1000.bp..42607752.5.Total..normalization.factor...0.23..effective.total...10000000.")]
MCF7_anno_Enhancer_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> MCF7_anno_Enhancer_strengh
colnames(MCF7_anno_Enhancer_strengh) <- c("symbol","Count_enh")

MCF7_anno_promoter_strengh <- MCF7_anno_Promoter[,c("Gene.Name","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_MC7..Tag.Count.in.1000.bp..42607752.5.Total..normalization.factor...0.23..effective.total...10000000.")]
MCF7_anno_promoter_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> MCF7_anno_promoter_strengh
colnames(MCF7_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(MCF7_anno_Enhancer_strengh, MCF7_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh

# Add gene expression information to the genes with nearby MYC binding if the gene is significant regulated upon MYCMI6 treatment in MCF7. The dataset is downloaded from GSE197062ID
plot2 <- read.table("DE_analysis_MYCMI6_MCF7.txt")
plot2 <- plot2[plot2$padj.c1 < .05,]
plot2 <- plot2[,c("symbol", "log2FC.c1")]
plot3 <- left_join(plot, plot2)
plot3[is.na(plot3)] = 0
plot3 <- plot3[!duplicated(plot3$symbol),]

# Load in ABC analysis which has been performed with HI-C from GSE66733
ABC_MCF7 <- read.delim("MCF7_MED_EnhancerPredictionsFull.txt")
ABC_MCF7_MYC <- read.delim("EnhancerPredictionsFull_MYC_binding_MCF7.txt")

ABC_merged <- merge(ABC_MCF7, ABC_MCF7_MYC, by.x = "name", by.y = "PeakID..cmd.annotatePeaks.pl.MCF7_Predictions_MED1.MCF7_MED_EnhancerPredictionsFull.txt.hg38..d....CellLine_MYC_MED1_ChIP.tagdir.Combined_MYC_MC7...cpu.10.")

ABC_merged <- ABC_merged[ABC_merged$...CellLine_MYC_MED1_ChIP.tagdir.Combined_MYC_MC7..Tag.Count.in.given.bp..42607752.5.Total..normalization.factor...0.23..effective.total...10000000. > 30,]

ABC_plot <- ABC_merged[,c("TargetGene","ABC.Score","powerlaw_contact")]
ABC_plot %>% group_by(TargetGene) %>% summarise_all(funs(sum)) -> ABC_plot_1
colnames(ABC_plot_1) <- c("symbol","ABC.Score","HIC")


plot4 <- left_join(plot3, ABC_plot_1)
plot4[is.na(plot4)] = 0

# Summarize bins for plotting
plot5 <- plot4[1:16350,2:7]
plot5 <- plot5[order(plot5$count),]
plot6 <- rowsum(plot5,rep(1:327,each=50))

ggplot() +
  geom_smooth(data= plot6 ,mapping=aes(x=rank(-count),y=Count_enh/50)) +
  geom_smooth(data= plot6 ,mapping=aes(x=rank(-count),y=Count_prom/50)) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + xlab("Ranked Promoter - Enhancer dominace") + ylab("MYC tag density") +
  theme(plot.title = element_text(hjust = 0.5))

# Density for ABC sscore over enhancer promoter ranking
ggplot() +
  geom_smooth(data= plot6 ,mapping=aes(x=rank(-count),y=(ABC.Score/50))) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + xlab("Ranked genes") + ylab("ABC score") +
  theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(ylim=c(0,0.08))

################# Supp 2E ################# 
MCF7_count <- read.delim("MYC_MCF7_Consensus.count.txt")

# Identify MYC promoters and enhancers
MCF7_count$Ann <- sub(" .*","",MCF7_count$Annotation)
MCF7_anno_Enhancer <- MCF7_count[MCF7_count$Ann == "Intergenic" | MCF7_count$Ann == "intron",]
MCF7_anno_Enhancer <- MCF7_anno_Enhancer[(MCF7_anno_Enhancer$Distance.to.TSS < 50000 & MCF7_anno_Enhancer$Detailed.Annotation > 3000) & MCF7_anno_Enhancer$Distance.to.TSS > -50000,]
MCF7_anno_Promoter <- MCF7_count[MCF7_count$Ann == "promoter-TSS",]
MCF7_anno_Promoter <- MCF7_anno_Promoter[!duplicated(MCF7_anno_Promoter$Gene.Name),]

# Summarise enhancer contribution 
MCF7_anno_Enhancer_strengh <- MCF7_anno_Enhancer[,c("Gene.Name","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_MC7..Tag.Count.in.1000.bp..42607752.5.Total..normalization.factor...0.23..effective.total...10000000.")]
MCF7_anno_Enhancer_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> MCF7_anno_Enhancer_strengh
colnames(MCF7_anno_Enhancer_strengh) <- c("symbol","Count_enh")

# Summarise promoter contribution 
MCF7_anno_promoter_strengh <- MCF7_anno_Promoter[,c("Gene.Name","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_MC7..Tag.Count.in.1000.bp..42607752.5.Total..normalization.factor...0.23..effective.total...10000000.")]
MCF7_anno_promoter_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> MCF7_anno_promoter_strengh
colnames(MCF7_anno_promoter_strengh) <- c("symbol","Count_prom")

# Merge genes with either promoter or enhancer MYC contribution
plot <- full_join(MCF7_anno_Enhancer_strengh, MCF7_anno_promoter_strengh)
plot[is.na(plot)] = 0

# Calculate the difference en contribution from different regions
plot$count <- plot$Count_prom-plot$Count_enh

# Add gene expression information to the genes with nearby MYC binding if the gene is significant regulated upon siRNA knockdown (data from DESeq2 script)
plot2 <- read.table("DE_analysis_MYCMI6_MCF7.txt")
plot2 <- plot2[plot2$padj.c1 < .05,]
plot2 <- plot2[,c("symbol", "log2FC.c1")]
plot3 <- left_join(plot, plot2)
plot3[is.na(plot3)] = 0
plot3 <- plot3[!duplicated(plot3$symbol),]

plot_sig <- plot3[abs(plot3$log2FC.c1) > 0,]

plot_prom <- plot_sig[plot_sig$count>quantile(plot_sig$count, probs=0.75),"symbol"]
plot_enh <- plot_sig[-plot_sig$count>quantile(-plot_sig$count, probs=0.75),"symbol"]

plot_prom <- as.data.frame(mapIds(org.Hs.eg.db, plot_prom$symbol, "ENTREZID",'SYMBOL'))
plot_enh <- as.data.frame(mapIds(org.Hs.eg.db, plot_enh$symbol, "ENTREZID",'SYMBOL'))

MCF7_prom <- enricher(plot_prom$`mapIds(org.Hs.eg.db, plot_prom$symbol, "ENTREZID", "SYMBOL")`, TERM2GENE=m_t2g)
MCF7_prom <- as.data.frame(MCF7_prom@result)
MCF7_prom_sig <- MCF7_prom[MCF7_prom$p.adjust < 0.05,]

ggplot(data=MCF7_prom_sig, aes(x=reorder(ID, -log10(p.adjust)), y=-log10(p.adjust))) +geom_bar(stat="identity", position = "dodge") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + coord_flip() + xlab(" ") + ylab("-log10(p.adjust)")


MCF7_enh <- enricher(plot_enh$`mapIds(org.Hs.eg.db, plot_enh$symbol, "ENTREZID", "SYMBOL")`, TERM2GENE=m_t2g)
MCF7_enh <- as.data.frame(MCF7_enh@result)
MCF7_enh_sig <- MCF7_enh[MCF7_enh$p.adjust < 0.05,]

ggplot(data=MCF7_enh_sig, aes(x=reorder(ID, -log10(p.adjust)), y=-log10(p.adjust))) +geom_bar(stat="identity", position = "dodge") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + scale_fill_viridis_c(option="E") + coord_flip() + xlab(" ") + ylab("-log10(p.adjust)")

################# Supp 2F #################  
MCF7_count <- read.delim("ER_MCF7_ann_ER_count.txt")

# Identify MYC promoters and enhancers
MCF7_count$Ann <- sub(" .*","",MCF7_count$Annotation)
MCF7_anno_Enhancer <- MCF7_count[MCF7_count$Ann == "Intergenic" | MCF7_count$Ann == "intron",]
MCF7_anno_Enhancer <- MCF7_anno_Enhancer[(MCF7_anno_Enhancer$Distance.to.TSS < 50000 & MCF7_anno_Enhancer$Detailed.Annotation > 3000) & MCF7_anno_Enhancer$Distance.to.TSS > -50000,]
MCF7_anno_Promoter <- MCF7_count[MCF7_count$Ann == "promoter-TSS",]
MCF7_anno_Promoter <- MCF7_anno_Promoter[!duplicated(MCF7_anno_Promoter$Gene.Name),]

# Summarise enhancer contribution 
MCF7_anno_Enhancer_strengh <- MCF7_anno_Enhancer[,c("Gene.Name","ER_MCF7.sorted.dup.bam.TD.Tag.Count.in.1000.bp..8265958.0.Total..normalization.factor...1.21..effective.total...10000000.")]
MCF7_anno_Enhancer_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> MCF7_anno_Enhancer_strengh
colnames(MCF7_anno_Enhancer_strengh) <- c("symbol","Count_enh")

# Summarise promoter contribution 
MCF7_anno_promoter_strengh <- MCF7_anno_Promoter[,c("Gene.Name","ER_MCF7.sorted.dup.bam.TD.Tag.Count.in.1000.bp..8265958.0.Total..normalization.factor...1.21..effective.total...10000000.")]
MCF7_anno_promoter_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> MCF7_anno_promoter_strengh
colnames(MCF7_anno_promoter_strengh) <- c("symbol","Count_prom")

# Merge genes with either promoter or enhancer MYC contribution
plot <- full_join(MCF7_anno_Enhancer_strengh, MCF7_anno_promoter_strengh)
plot[is.na(plot)] = 0

# Calculate the difference en contribution from different regions
plot$count <- plot$Count_prom-plot$Count_enh

# Add gene expression information to the genes with nearby MYC binding if the gene is significant regulated upon siRNA knockdown (data from DESeq2 script)
plot2 <- read.table("DE_analysis_KO_MCF7.txt")
plot2 <- plot2[plot2$padj_MYC_vs_Crtl < .05,]
plot2 <- plot2[,c("symbol", "Log2FC_MYC_vs_Crtl")]
plot3 <- left_join(plot, plot2)
plot3[is.na(plot3)] = 0
plot3 <- plot3[!duplicated(plot3$symbol),]

# Summarize bins for plotting
plot5 <- plot3[1:6200,2:5]
plot5 <- plot5[order(plot5$count),]
plot6 <- rowsum(plot5,rep(1:124,each=50))

ggplot() +
  geom_smooth(data= plot6 ,mapping=aes(x=rank(-count),y=Count_enh/50), col = "red") +
  geom_smooth(data= plot6 ,mapping=aes(x=rank(-count),y=Count_prom/50), col = "blue") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + xlab("Ranked Promoter - Enhancer dominace") + ylab("ER tag density") +
  theme(plot.title = element_text(hjust = 0.5))

################# Supp 2G ################# 
MB231_count <- read.delim("STAT3_MB231_count.txt")

# Identify MYC promoters and enhancers
MB231_count$Ann <- sub(" .*","",MB231_count$Annotation)
MB231_anno_Enhancer <- MB231_count[MB231_count$Ann == "Intergenic" | MB231_count$Ann == "intron",]
MB231_anno_Enhancer <- MB231_anno_Enhancer[(MB231_anno_Enhancer$Distance.to.TSS < 50000 & MB231_anno_Enhancer$Detailed.Annotation > 3000) & MB231_anno_Enhancer$Distance.to.TSS > -50000,]
MB231_anno_Promoter <- MB231_count[MB231_count$Ann == "promoter-TSS",]
MB231_anno_Promoter <- MB231_anno_Promoter[!duplicated(MB231_anno_Promoter$Gene.Name),]

# Summarise enhancer contribution 
MB231_anno_Enhancer_strengh <- MB231_anno_Enhancer[,c("Gene.Name","STAT3_MDA231.sorted.dup.bam.TD.Tag.Count.in.given.bp..16476809.0.Total..normalization.factor...0.61..effective.total...10000000.")]
MB231_anno_Enhancer_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> MB231_anno_Enhancer_strengh
colnames(MB231_anno_Enhancer_strengh) <- c("symbol","Count_enh")

# Summarise promoter contribution 
MB231_anno_promoter_strengh <- MB231_anno_Promoter[,c("Gene.Name","STAT3_MDA231.sorted.dup.bam.TD.Tag.Count.in.given.bp..16476809.0.Total..normalization.factor...0.61..effective.total...10000000.")]
MB231_anno_promoter_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> MB231_anno_promoter_strengh
colnames(MB231_anno_promoter_strengh) <- c("symbol","Count_prom")

# Merge genes with either promoter or enhancer MYC contribution
plot <- full_join(MB231_anno_Enhancer_strengh, MB231_anno_promoter_strengh)
plot[is.na(plot)] = 0

# Calculate the difference en contribution from different regions
plot$count <- plot$Count_prom-plot$Count_enh

# Add gene expression information to the genes with nearby MYC binding if the gene is significant regulated upon siRNA knockdown (data from DESeq2 script)
plot2 <- read.table("DE_analysis_KO_MB231.txt")
plot2 <- plot2[plot2$padj_MYC_vs_Crtl < .05,]
plot2 <- plot2[,c("symbol", "Log2FC_MYC_vs_Crtl")]
plot3 <- left_join(plot, plot2)
plot3[is.na(plot3)] = 0
plot3 <- plot3[!duplicated(plot3$symbol),]

# Summarize bins for plotting
plot5 <- plot3[1:13250,2:5]
plot5 <- plot5[order(plot5$count),]
plot6 <- rowsum(plot5,rep(1:265,each=50))

ggplot() +
  geom_smooth(data= plot6 ,mapping=aes(x=rank(-count),y=Count_enh/50), col = "red") +
  geom_smooth(data= plot6 ,mapping=aes(x=rank(-count),y=Count_prom/50), col = "blue") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + xlab("Ranked Promoter - Enhancer dominace") + ylab("STAT3 tag density") +
  theme(plot.title = element_text(hjust = 0.5))


