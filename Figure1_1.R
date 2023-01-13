##############################################
################## FIGURE 1 ##################
##############################################

################# Figure 1A ################## 
# Annotation of MYC ChIP-seq binding in BT549 (consensus peaks from N=3)

# annotatePeaks.pl MYC_BT549_Consensus.bed hg38 > MYC_BT549_Consensus_ann.txt
annotation_BT549 <- read.delim("MYC_BT549_Consensus_ann.txt") # Load in dataframe
annotation_BT549$Ann <- sub(" .*","",annotation_BT549$Annotation) # Remove intron annotation
annotation_BT549_plot <- as.data.frame(table(annotation_BT549$Ann)) # Summaries data for plotting
annotation_BT549_plot$type <- rep("BT549")

annotation_MB231 <- read.delim("MYC_MB231_Consensus_ann.txt") # Load in dataframe
annotation_MB231$Ann <- sub(" .*","",annotation_MB231$Annotation) # Remove intron annotation
annotation_MB231_plot <- as.data.frame(table(annotation_MB231$Ann)) # Summaries data for plotting
annotation_MB231_plot$type <- rep("MB231")

annotation_MCF7 <- read.delim("MYC_MCF7_Consensus_ann.txt") # Load in dataframe
annotation_MCF7$Ann <- sub(" .*","",annotation_MCF7$Annotation) # Remove intron annotation
annotation_MCF7_plot <- as.data.frame(table(annotation_MCF7$Ann)) # Summaries data for plotting
annotation_MCF7_plot$type <- rep("MCF7")

annotation_T47D <- read.delim("MYC_T47D_Consensus_ann.txt") # Load in dataframe
annotation_T47D$Ann <- sub(" .*","",annotation_T47D$Annotation) # Remove intron annotation
annotation_T47D_plot <- as.data.frame(table(annotation_T47D$Ann)) # Summaries data for plotting
annotation_T47D_plot$type <- rep("T47D")

plot <- rbind(annotation_BT549_plot, annotation_MB231_plot, annotation_MCF7_plot, annotation_T47D_plot)

library(tidyverse)
library(viridis)

ggplot(plot, aes(x=type, y=Freq, fill  = Var1)) +  geom_bar(position="fill", stat="identity") + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "right")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("% Binding distribution") + xlab(" ")

################# Figure 1B ################## 
annotation_BT549_int <- annotation_BT549[annotation_BT549$Ann == "Intergenic" |annotation_BT549$Ann == "intron",2:4]
annotation_MB231 <- annotation_MB231[annotation_MB231$Ann == "Intergenic" |annotation_MB231$Ann == "intron",2:4]
annotation_MCF7 <- annotation_MCF7[annotation_MCF7$Ann == "Intergenic" |annotation_MCF7$Ann == "intron",2:4]
annotation_T47D <- annotation_T47D[annotation_T47D$Ann == "Intergenic" |annotation_T47D$Ann == "intron",2:4]

enhancer <- rbind(annotation_BT549_int,annotation_MB231,annotation_T47D,annotation_MCF7)

write.table(enhancer, "BC_MYC_Commonenhancer.bed", sep = "\t", quote = F, row.names = F)

data <- c(53805,12036)
info <- c("Encode","Unknown")

plot <- as.data.frame(rbind(data,info))

ggplot(plot, aes(x="", y=data, fill=info)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "right")

################# Figure 1C ################## 
# annotatePeaks.pl MYC_Consensus.bed hg38 -d tagdir/Combined_MED1_BT549 tagdir/Combined_MED1_MCF7 tagdir/Combined_MED1_MDA231 tagdir/Combined_MED1_MDA468 tagdir/Combined_MED1_T47D tagdir/Combined_MYC_BT549 tagdir/Combined_MYC_MCF7 tagdir/Combined_MYC_MDA231 tagdir/Combined_MYC_MDA468 tagdir/Combined_MYC_T47D
library(corrplot)
MYC_Counts <- read.delim("MYC_Concensus.count.txt") # Load matrix of all identified MYC binding sites
MYC_Counts$Ann <- sub(" .*","",MYC_Counts$Annotation) # Change annnotation for subsetting
MYC_Counts_Intergenic <- MYC_Counts[MYC_Counts$Ann == "Intergenic" | MYC_Counts$Ann == "intron",] # Get only intergenic regions
MYC_Counts_Intergenic <- MYC_Counts_Intergenic[,20:29] # subset counts for each data set
colnames(MYC_Counts_Intergenic) <- c("MED1_BT549","MED1_MCF7","MED1_MDA231","MED1_MDA468","MED1_T47D","MYC_BT549","MYC_MCF7","MYC_MDA231","MYC_MDA468","MYC_T47D") # Change colnames to something easy understandable
t1 <- MYC_Counts_Intergenic[,c("MYC_BT549","MYC_MDA231","MYC_T47D","MYC_MCF7")] # Prepare for heatmap plot
res_intergenic <- cor(t1, method = c("spearman")) # Calculate correlation
corrplot(res_intergenic, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45,  col.lim = c(0, 1), is.corr = FALSE, col=colorRampPalette(c("white","#B11518"))(200) ,tl.cex = 2) # plot correlation matrix

MYC_Counts_Promoter <- MYC_Counts[MYC_Counts$Ann == "promoter-TSS",]
MYC_Counts_Promoter <- MYC_Counts_Promoter[,20:29]
colnames(MYC_Counts_Promoter) <- c("MED1_BT549","MED1_MCF7","MED1_MDA231","MED1_MDA468","MED1_T47D","MYC_BT549","MYC_MCF7","MYC_MDA231","MYC_MDA468","MYC_T47D")
t1 <- MYC_Counts_Promoter[,c("MYC_BT549","MYC_MDA231","MYC_T47D","MYC_MCF7")]
res_promoter <- cor(t1, method = c("spearman"))
corrplot(res_promoter, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45,  col.lim = c(0, 1), is.corr = FALSE, col=colorRampPalette(c("white","#253884"))(200) ,tl.cex = 2) 


################# Figrue 1D ################# 
library(corrplot)
MYC_Counts <- read.delim("MYC_Master_ann.txt") # Load matrix of all identified MYC binding sites
MYC_Counts$Ann <- sub(" .*","",MYC_Counts$Annotation) # Change annnotation for subsetting
MYC_Counts_Intergenic <- MYC_Counts[MYC_Counts$Ann == "Intergenic",] # Get only intergenic regions
MYC_Counts_Intergenic <- MYC_Counts_Intergenic[,20:32] # subset counts for each data set
colnames(MYC_Counts_Intergenic) <- c("BT549","MCF7","A493","BT186","H128","H2171", "HCT116","HeLA","HT1080","MCF7.1", "MM1S", "U2OS", "U87") # Change colnames to something easy understandable
t1 <- MYC_Counts_Intergenic[,c("BT549","MCF7","A493","BT186","H128","H2171", "HCT116","HeLA","HT1080", "MM1S", "U2OS", "U87")] # Prepare for heatmap plot
res_intergenic <- cor(t1, method = c("spearman")) # Calculate correlation
corrplot(res_intergenic, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45,  col.lim = c(0, 1), is.corr = FALSE, col=colorRampPalette(c("white","#B11518"))(200),tl.cex = 2) # plot correlation matrix

MYC_Counts_Promoter <- MYC_Counts[MYC_Counts$Ann == "promoter-TSS",]
MYC_Counts_Promoter <- MYC_Counts_Promoter[,20:32]
colnames(MYC_Counts_Promoter) <- c("BT549","MCF7","A493","BT186","H128","H2171", "HCT116","HeLA","HT1080","MCF7.1", "MM1S", "U2OS", "U87")
t1 <- MYC_Counts_Promoter[,c("BT549","MCF7","A493","BT186","H128","H2171", "HCT116","HeLA","HT1080", "MM1S", "U2OS", "U87")]
res_promoter <- cor(t1, method = c("spearman"))
corrplot(res_promoter, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45,  col.lim = c(0, 1), is.corr = FALSE,  col=colorRampPalette(c("white","#253884"))(200),tl.cex = 2) 




################# Figure 1E ################# 

# DESeq2 analysis of enhancers
library(DESeq2)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(pheatmap)
library(viridis)

##Differential gene expression##
#Import count table generated by featureCounts and make metadata
a <- read.delim("MYC_Consensus.raw.counts.txt")
a2<-a[,39:53] #Removes chr, start, end for exons which take up a lot of space
colnames(a2) <- c("BT549_1","BT549_2","BT549_3","MCF7_1","MCF7_2","MCF7_3","MDA231_1","MDA231_2","MDA231_3","MDA468_1","MDA468_2","MDA468_3","T47D_1","T47D_2","T47D_3")
a2 <- a2[,c("BT549_1","BT549_2","BT549_3","MCF7_1","MCF7_2","MCF7_3","MDA231_1","MDA231_2","MDA231_3","T47D_2","T47D_3")]

# Make metadata
grp1 <- c("BT549_1","BT549_2","BT549_3")
grp2 <- c("MCF7_1","MCF7_2","MCF7_3")
grp3 <- c("MDA231_1","MDA231_2","MDA231_3")
grp4 <- c("T47D_2","T47D_3")
cond1="B5549"
cond2="MCF7"
cond3="MDA231"
cond4="T47D"
condition=factor(c(rep(cond1,length(grp1)),rep(cond2,length(grp2)),rep(cond3,length(grp3)),rep(cond4,length(grp4))))
condition <- relevel(condition, ref=cond1)
replicate<-factor(c("a","b","c","a","b","c","a","b","c","b","c"))

# Mate design tabel for potential batch normalization
sampleTable <- data.frame(
  row.names=c(grp1,grp2,grp3,grp4),
  replicate=replicate,
  condition=condition)
sampleTable

n=6 #Set number of comparisons. This will be used for deseq2. The list below determines the actual comparisons. 
conList<-list(c("condition","B5549","MCF7"),
              c("condition","B5549","MDA231"),
              c("condition","B5549","T47D"),
              c("condition","MCF7","MDA231"),
              c("condition","MCF7","T47D"),
              c("condition","MDA231","T47D"))

# Annotate with gene names
gene <- a$PeakID..cmd.annotatePeaks.pl.peak.MYC_Consensus.bed.hg38..d.tagdir.BT549_input_pool_9094_S27_.sorted.dup.bam.TD.tagdir.MCF7_input_pool_9091_S24_.sorted.dup.bam.TD.tagdir.MDA231_input_pool_9095_S28_.sorted.dup.bam.TD.tagdir.MDA468_input_pool_9093_S26_.sorted.dup.bam.TD.tagdir.MED1_BT549_R1_9098_S31_.sorted.dup.bam.TD.tagdir.MED1_BT549_R2_9099_S32_.sorted.dup.bam.TD.tagdir.MED1_BT549_R3_9084_S17_.sorted.dup.bam.TD.tagdir.MED1_MCF7_R1_9077_S10_.sorted.dup.bam.TD.tagdir.MED1_MCF7_R2_9096_S29_.sorted.dup.bam.TD.tagdir.MED1_MCF7_R3_9097_S30_.sorted.dup.bam.TD.tagdir.MED1_MDA231_R1_9085_S18_.sorted.dup.bam.TD.tagdir.MED1_MDA231_R2_9086_S19_.sorted.dup.bam.TD.tagdir.MED1_MDA231_R3_9100_S33_.sorted.dup.bam.TD.tagdir.MED1_MDA468_R1_9081_S14_.sorted.dup.bam.TD.tagdir.MED1_MDA468_R2_9082_S15_.sorted.dup.bam.TD.tagdir.MED1_MDA468_R3_9083_S16_.sorted.dup.bam.TD.tagdir.MED1_T47D_R1_9078_S11_.sorted.dup.bam.TD.tagdir.MED1_T47D_R2_9079_S12_.sorted.dup.bam.TD.tagdir.MED1_T47D_R3_9080_S13_.sorted.dup.bam.TD.tagdir.MYC_BT549_R1_9104_S37_.sorted.dup.bam.TD.tagdir.MYC_BT549_R2_9105_S38_.sorted.dup.bam.TD.tagdir.MYC_BT549_R3_9106_S39_.sorted.dup.bam.TD.tagdir.MYC_MCF7_R1_9101_S34_.sorted.dup.bam.TD.tagdir.MYC_MCF7_R2_9108_S41_.sorted.dup.bam.TD.tagdir.MYC_MCF7_R3_9109_S42_.sorted.dup.bam.TD.tagdir.MYC_MDA231_R1_9111_S44_.sorted.dup.bam.TD.tagdir.MYC_MDA231_R2_9107_S40_.sorted.dup.bam.TD.tagdir.MYC_MDA231_R3_9090_S23_.sorted.dup.bam.TD.tagdir.MYC_MDA468_R1_9103_S36_.sorted.dup.bam.TD.tagdir.MYC_MDA468_R2_9088_S21_.sorted.dup.bam.TD.tagdir.MYC_MDA468_R3_9089_S22_.sorted.dup.bam.TD.tagdir.MYC_T47D_R1_9087_S20_.sorted.dup.bam.TD.tagdir.MYC_T47D_R2_9102_S35_.sorted.dup.bam.TD.tagdir.MYC_T47D_R3_9110_S43_.sorted.dup.bam.TD.tagdir.T47D_input_pool_9092_S25_.sorted.dup.bam.TD..size.1000..raw.

#Make differential analysis using DESeq2
library(DESeq2)
a_final <- DESeqDataSetFromMatrix(countData = round(a2), colData = sampleTable, design = ~condition)
a_final <- DESeq(a_final) #Run diffexp analysis


#Attach padj and logFC for each comparison to raw data
a6 <- gene_count <- as.data.frame(counts(a_final, normalized=TRUE))
a6 <- cbind(a6, a[,2:19])
diff<-a6 # Normalized counts
for (i in 1:n) {
  res<-data.frame(results(a_final,contrast = conList[[i]])) #Get p values for different comparisons 
  colnames(res)[c(2,6)]<-c(paste("log2FC",paste(".c",i,sep = ""),sep = ""),paste("padj",paste(".c",i,sep = ""),sep = ""))
  diff<-cbind(diff,res[,c(2,6)])
}

# Make row names nice!
colnames(diff) <- c("BT549_1","BT549_2","BT549_3","MCF7_1","MCF7_2","MCF7_3","MDA231_1","MDA231_2","MDA231_3","T47D_2","T47D_3","Chr","Start","End","Strand","Peak.Score","Region.Size","Annotation","Detailed.Annotation","Distances.to.TSS","Nearest.PromoterID","Entrez.ID","Nearest.Unigene","Nearest.Refseq","Nearest.Ensembl","Gene.Name","Gene.Alias","Gene.Descriotion","Gene.Type",
                    "BT549_vs_MCF7_log","BT549_vs_MCF7_padj","BT549_vs_MDA231_log","BT549_vs_MDA231_padj","BT549_vs_T47D_log","BT549_vs_T47D_padj","MCF7_vs_MDA231_log","MCF7_vs_MDA231_padj","MCF7_vs_T47D_log","MCF7_vs_T47D_padj","MDA231_vs_T47D_log","MDA231_vs_T47D_padj")


# Devide into based on promoters and enhancers 
diff$Ann <- sub(" .*","",diff$Annotation)
diff_intergenic <- diff[diff$Ann == "Intergenic" | diff$Ann == "intron",]
diff_promoter <- diff[diff$Ann == "promoter-TSS",]

# Define MESO specific enhancer
diff_intergenic_Meso <- diff_intergenic[diff_intergenic$BT549_vs_MCF7_log > 0 & diff_intergenic$BT549_vs_T47D_log > 0 & diff_intergenic$MCF7_vs_MDA231_log < 0 & diff_intergenic$MDA231_vs_T47D_log > 0 & 
                                          diff_intergenic$BT549_vs_MCF7_padj < 0.1 & diff_intergenic$BT549_vs_T47D_padj < 0.1 & diff_intergenic$MCF7_vs_MDA231_padj < 0.1 & diff_intergenic$MDA231_vs_T47D_padj < 0.1,]
diff_intergenic_Recep <- diff_intergenic[diff_intergenic$BT549_vs_MCF7_log < 0 & diff_intergenic$BT549_vs_T47D_log < 0 & diff_intergenic$MCF7_vs_MDA231_log > 0 & diff_intergenic$MDA231_vs_T47D_log < 0 & 
                                           diff_intergenic$BT549_vs_MCF7_padj < 0.1 & diff_intergenic$BT549_vs_T47D_padj < 0.1 & diff_intergenic$MCF7_vs_MDA231_padj < 0.1 & diff_intergenic$MDA231_vs_T47D_padj < 0.1,]

boxplot(diff_intergenic_Meso$BT549_1, diff_intergenic_Meso$MDA231_1, diff_intergenic_Meso$MCF7_1, diff_intergenic_Meso$T47D_2, outline = F)
boxplot(diff_intergenic_Recep$BT549_1, diff_intergenic_Recep$MDA231_1, diff_intergenic_Recep$MCF7_1, diff_intergenic_Recep$T47D_2, outline = F)

write.table(diff_intergenic_Meso[12:14], "Meso_MYC_Enh.nogene.bed", sep = "\t", quote = F, row.names = F)
write.table(diff_intergenic_Recep[12:14], "ER_MYC_Enh.nogene.bed", sep = "\t", quote = F, row.names = F)

# Make data with homer histogram and plot
library(pheatmap)
library(viridis)
library(corrplot)
library(ComplexHeatmap)
library(circlize)
library(magick)
library(reshape2)
library(ggsci)
library(VennDiagram)
a<-read.delim("Meso_MYC_Enh_ghist1.txt",h=T) 
a2<-a[order(a$X0.5+a$X0.7, decreasing = T),2:ncol(a)] 
MED1_BT549 <-a2[,1:101]
MED1_MCF7 <-a2[,102:202]
MED1_MDA231 <-a2[,203:303]
MED1_MDA468 <-a2[,304:404]
MED1_T47D <-a2[,405:505]
MYC_BT549 <-a2[,506:606]
MYC_MCF7 <-a2[,607:707]
MYC_MDA231 <-a2[,708:808]
MYC_MDA468 <-a2[,809:909]
MYC_T47D <-a2[,910:1010]
ER_MCF7 <-a2[,1011:1111]
STAT3_MDA231 <-a2[,1112:1212]


BT549_MED <- Heatmap(as.matrix(MED1_BT549), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("BT549", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MCF7_MED <- Heatmap(as.matrix(MED1_MCF7), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MCF7", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MDA231_MED <- Heatmap(as.matrix(MED1_MDA231), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MDA231", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(12,"cm"), width = unit(1.5,"cm"))
T47D_MED <- Heatmap(as.matrix(MED1_T47D), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("T47D", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
BT549_MYC <- Heatmap(as.matrix(MYC_BT549), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","green1","green4")), use_raster = T, raster_quality = 10, column_split = c(rep("BT549", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MCF7_MYC <- Heatmap(as.matrix(MYC_MCF7), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","green1","green4")), use_raster = T, raster_quality = 10, column_split = c(rep("MCF7", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MDA231_MYC <- Heatmap(as.matrix(MYC_MDA231), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","green1","green4")), use_raster = T, raster_quality = 10, column_split = c(rep("MDA231", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
T47D_MYC <- Heatmap(as.matrix(MYC_T47D), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","green1","green4")), use_raster = T, raster_quality = 10, column_split = c(rep("T47D", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MCF7 <- Heatmap(as.matrix(ER_MCF7), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","red1","red4")), use_raster = T, raster_quality = 10, column_split = c(rep("ER_MCF7", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MDA231 <- Heatmap(as.matrix(STAT3_MDA231), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","red1","red4")), use_raster = T, raster_quality = 10, column_split = c(rep("STAT3_MDA231", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))

draw(BT549_MYC+MDA231_MYC+MCF7_MYC+T47D_MYC+BT549_MED+MDA231_MED+MCF7_MED+T47D_MED+MCF7+MDA231, heatmap_legend_side="bottom",row_title=paste(nrow(a2),"sites", collapse = " "))



library(pheatmap)
library(viridis)
library(corrplot)
library(ComplexHeatmap)
library(circlize)
library(magick)
library(reshape2)
library(ggsci)
library(VennDiagram)
a<-read.delim("ER_MYC_Enh_ghist1.txt",h=T) 
a2<-a[order(a$X0.6+a$X0.8, decreasing = T),2:ncol(a)] 
MED1_BT549 <-a2[,1:101]
MED1_MCF7 <-a2[,102:202]
MED1_MDA231 <-a2[,203:303]
MED1_MDA468 <-a2[,304:404]
MED1_T47D <-a2[,405:505]
MYC_BT549 <-a2[,506:606]
MYC_MCF7 <-a2[,607:707]
MYC_MDA231 <-a2[,708:808]
MYC_MDA468 <-a2[,809:909]
MYC_T47D <-a2[,910:1010]
ER_MCF7 <-a2[,1011:1111]
STAT3_MDA231 <-a2[,1112:1212]


BT549_MED <- Heatmap(as.matrix(MED1_BT549), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("BT549", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MCF7_MED <- Heatmap(as.matrix(MED1_MCF7), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MCF7", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MDA231_MED <- Heatmap(as.matrix(MED1_MDA231), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MDA231", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(12,"cm"), width = unit(1.5,"cm"))
T47D_MED <- Heatmap(as.matrix(MED1_T47D), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("T47D", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
BT549_MYC <- Heatmap(as.matrix(MYC_BT549), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","green1","green4")), use_raster = T, raster_quality = 10, column_split = c(rep("BT549", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MCF7_MYC <- Heatmap(as.matrix(MYC_MCF7), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","green1","green4")), use_raster = T, raster_quality = 10, column_split = c(rep("MCF7", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MDA231_MYC <- Heatmap(as.matrix(MYC_MDA231), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","green1","green4")), use_raster = T, raster_quality = 10, column_split = c(rep("MDA231", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
T47D_MYC <- Heatmap(as.matrix(MYC_T47D), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","green1","green4")), use_raster = T, raster_quality = 10, column_split = c(rep("T47D", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MCF7 <- Heatmap(as.matrix(ER_MCF7), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,5,10),c("white","red1","red4")), use_raster = T, raster_quality = 10, column_split = c(rep("ER_MCF7", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MDA231 <- Heatmap(as.matrix(STAT3_MDA231), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,5,10),c("white","red1","red4")), use_raster = T, raster_quality = 10, column_split = c(rep("STAT3_MDA231", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))

draw(BT549_MYC+MDA231_MYC+MCF7_MYC+T47D_MYC+BT549_MED+MDA231_MED+MCF7_MED+T47D_MED+MCF7+MDA231, heatmap_legend_side="bottom",row_title=paste(nrow(a2),"sites", collapse = " "))


################# Figure 1F ################# 

# mergePeaks /data/Siersbaek/STJ_RS/MED1_lost_sites_RM46_Veh_Tagthreshold2.pos BRCA_log2norm.bed > RM_MED1_Lost.txt
# grep "/data/Siersbaek/STJ_RS/MED1_lost_sites_RM46_Veh_Tagthreshold2.pos|BRCA_log2norm.bed" RM_MED1_Lost.txt > RM_MED1_Lost.bed
# awk '{print $10}' RM_MED1_Lost.bed > RM_MED1_Lost.list
# grep -w -f RM_MED1_Lost.list BRCA_log2norm.txt > RM_MED1_Lost.Count.txt
# cat BRCA_names.txt RM_MED1_Lost.Count.txt > RM_MED1_Lost.Count1.txt
Meso_ATAC <- read.delim("Meso_MYC_BRCA_log2norm1.txt")
ER_ATAC <- read.delim("ER_MYC_BRCA_log2norm1.txt")

ATAC_Info <- read.delim("TCGA_identifier_mapping.txt")
ATAC_Info$identifier <- gsub("-", "_", ATAC_Info$bam_prefix)

ATAC_Into_BC <- merge(ATAC_Info, tcga_clin, by.x = "patient_id", by.y = "PATIENT_ID")

ATAC_TNBC <- ATAC_Into_BC[ATAC_Into_BC$ER_STATUS_BY_IHC == "Negative" & ATAC_Into_BC$PR_STATUS_BY_IHC == "Negative",]
ATAC_ER  <- ATAC_Into_BC[ATAC_Into_BC$ER_STATUS_BY_IHC == "Positive",]

ATAC_TNBC_count_TNBC <- Meso_ATAC[,colnames(Meso_ATAC) %in% ATAC_TNBC$identifier]
ATAC_TNBC_count_ER <- Meso_ATAC[,colnames(Meso_ATAC) %in% ATAC_ER$identifier]

ATAC_ER_count_TNBC <- ER_ATAC[,colnames(ER_ATAC) %in% ATAC_TNBC$identifier]
ATAC_ER_count_ER <- ER_ATAC[,colnames(ER_ATAC) %in% ATAC_ER$identifier]

Meso_TN <- as.data.frame(colMeans(ATAC_TNBC_count_TNBC))
Meso_ER <- as.data.frame(colMeans(ATAC_TNBC_count_ER))
ER_TN <- as.data.frame(colMeans(ATAC_ER_count_TNBC))
ER_ER <- as.data.frame(colMeans(ATAC_ER_count_ER))


Meso_TN$type <- rep("Meso_Basal")
Meso_ER$type <- rep("Meso_ER")
ER_TN$type <- rep("ER_Basal")
ER_ER$type <- rep("ER_ER")

colnames(Meso_TN) <- c("log","type")
colnames(Meso_ER) <- c("log","type")
colnames(ER_TN) <- c("log","type")
colnames(ER_ER) <- c("log","type")

plot <- rbind(Meso_TN, Meso_ER, ER_TN, ER_ER)
plot$type <- factor(plot$type, c("Meso_Basal","Meso_ER","ER_Basal","ER_ER"))

ggplot(plot, aes(x=type, y=log, fill = type)) + geom_violin(outlier.alpha = 1, alpha = 0.8) + geom_boxplot(alpha = 0.4, outlier.alpha = 0) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  xlab(" ") + ylab("Log(FPKM)") + scale_fill_viridis_d(option="E") + coord_cartesian(ylim=c(0.5,3))


### same for promoters
Meso_ATAC <- read.delim("Meso_Prom1.txt")
ER_ATAC <- read.delim("ER_Prom1.txt")

ATAC_Info <- read.delim("TCGA_identifier_mapping.txt")
ATAC_Info$identifier <- gsub("-", "_", ATAC_Info$bam_prefix)

ATAC_Into_BC <- merge(ATAC_Info, tcga_clin, by.x = "patient_id", by.y = "PATIENT_ID")

ATAC_TNBC <- ATAC_Into_BC[ATAC_Into_BC$ER_STATUS_BY_IHC == "Negative" & ATAC_Into_BC$PR_STATUS_BY_IHC == "Negative",]
ATAC_ER  <- ATAC_Into_BC[ATAC_Into_BC$ER_STATUS_BY_IHC == "Positive",]

ATAC_TNBC_count_TNBC <- Meso_ATAC[,colnames(Meso_ATAC) %in% ATAC_TNBC$identifier]
ATAC_TNBC_count_ER <- Meso_ATAC[,colnames(Meso_ATAC) %in% ATAC_ER$identifier]

ATAC_ER_count_TNBC <- ER_ATAC[,colnames(ER_ATAC) %in% ATAC_TNBC$identifier]
ATAC_ER_count_ER <- ER_ATAC[,colnames(ER_ATAC) %in% ATAC_ER$identifier]

Meso_TN <- as.data.frame(colMeans(ATAC_TNBC_count_TNBC))
Meso_ER <- as.data.frame(colMeans(ATAC_TNBC_count_ER))
ER_TN <- as.data.frame(colMeans(ATAC_ER_count_TNBC))
ER_ER <- as.data.frame(colMeans(ATAC_ER_count_ER))


Meso_TN$type <- rep("Meso_Basal")
Meso_ER$type <- rep("Meso_ER")
ER_TN$type <- rep("ER_Basal")
ER_ER$type <- rep("ER_ER")

colnames(Meso_TN) <- c("log","type")
colnames(Meso_ER) <- c("log","type")
colnames(ER_TN) <- c("log","type")
colnames(ER_ER) <- c("log","type")

plot <- rbind(Meso_TN, Meso_ER, ER_TN, ER_ER)
plot$type <- factor(plot$type, c("Meso_Basal","Meso_ER","ER_Basal","ER_ER"))

ggplot(plot, aes(x=type, y=log, fill = type)) + geom_violin(outlier.alpha = 1, alpha = 0.8) + geom_boxplot(alpha = 0.4, outlier.alpha = 0) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  xlab(" ") + ylab("Log(FPKM)") + scale_fill_viridis_d(option="E") + coord_cartesian(ylim=c(0.5,5))







### Supp 1C motif scores weak and strong
Enhancer <- read.delim("MYC_Consensus_annotated_Promoter_motif_hist.txt") 
avr2_enh <-Enhancer[,c(1,grep("total.sites", colnames(Enhancer)))]
colnames(avr2_enh)<-c("dist","Promoter_1","Promoter_2","Promoter")

Enhancer <- read.delim("MYC_Consensus_annotated_Intergenic_motif_hist.txt")
avr1_enh <-Enhancer[,c(1,grep("total.sites", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist_1","Enhancer_1","Enhancer_2","Enhancer")

plot <- cbind(avr1_enh,avr2_enh)
plot <- plot[,c("dist","Enhancer_1","Promoter_1")]

plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("Motif Enrichment") + xlab("bp") + xlim(-2000,2000) + ylim(0.0008,0.002)
