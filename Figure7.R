##############################################
################## FIGURE 7 ##################
##############################################

# Pre-proccessed data can be downloaded from https://zenodo.org/deposit/8323614

################# Packages #################
library(tidyverse)
library(PupillometryR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(reshape)
library(clusterProfiler)
library(pheatmap)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(magick)
library(ggsci)








################# Figure 7a #################
Enhancer <- read.delim("MED1_MYC_Enh_Down_GCN5_KDM3A_H3K9mod.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","GCN5_D","GCN5_K","KDM3A_D","KDM3A_HD","KDM3A_HK","KDM3A_K"
                      ,"Ac_D","Ac_HD","Ac_HK","Ac_K"
                      ,"Me1_D","Me1_HD","Me1_HK","Me1_K"
                      ,"Me2_D","Me2_HD","Me2_HK","Me2_K"
                      ,"Input_D","Input_HD","Input_HK","Input_K")

# GCN5 upon KJ in MYC activated enhancer
plot <- avr1_enh[,c("dist","GCN5_D","GCN5_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

# H3K9ac upon KJ in MYC activated enhancer
plot <- avr1_enh[,c("dist","Ac_D","Ac_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

# KDM3A upon KJ in MYC activated enhancer
plot <- avr1_enh[,c("dist","KDM3A_D","KDM3A_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

# H3K9me1 upon KJ in MYC activated enhancer
plot <- avr1_enh[,c("dist","Me1_D","Me1_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

# H3K9me2 upon KJ in MYC activated enhancer
plot <- avr1_enh[,c("dist","Me2_D","Me2_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

# annotatePeaks.pl tss hg38 -list ../BT549_KJ_POL/MYC_Gene_down.txt -d tagdig/GCN5_D.TD tagdig/GCN5_K.TD tagdig/KDM3A_D.TD tagdig/KDM3A_Hyp_D.TD tagdig/KDM3A_Hyp_K.TD tagdig/KDM3A_K.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_D.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/H3K9ac_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me1_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/H3K9me2_K.TD ../BT549_H3K9_KJ/tagdig/Input_D.TD ../BT549_H3K9_KJ/tagdig/Input_Hyp_D.TD ../BT549_H3K9_KJ/tagdig/Input_Hyp_K.TD ../BT549_H3K9_KJ/tagdig/Input_K.TD -size 5000 -hist 50 -cpu 12 > MED1_MYC_Gene_Down_GCN5_KDM3A_H3K9mod.txt
Promoter <- read.delim("MED1_MYC_Gene_Down_GCN5_KDM3A_H3K9mod.txt")
avr1_prom <-Promoter[,c(1,grep("Coverage", colnames(Promoter)))]
colnames(avr1_prom)<-c("dist","GCN5_D","GCN5_K","KDM3A_D","KDM3A_HD","KDM3A_HK","KDM3A_K"
                       ,"Ac_D","Ac_HD","Ac_HK","Ac_K"
                       ,"Me1_D","Me1_HD","Me1_HK","Me1_K"
                       ,"Me2_D","Me2_HD","Me2_HK","Me2_K"
                       ,"Input_D","Input_HD","Input_HK","Input_K")

# GCN5 upon KJ in MYC activated promoter/genes
plot <- avr1_prom[,c("dist","GCN5_D","GCN5_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

# GCN5 upon KJ in MYC activated promoter/genes
plot <- avr1_prom[,c("dist","Ac_D","Ac_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

# KDM3A upon KJ in MYC activated promoter/genes
plot <- avr1_prom[,c("dist","KDM3A_D","KDM3A_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

# H3K9me1 upon KJ in MYC activated promoter/genes
plot <- avr1_prom[,c("dist","Me1_D","Me1_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

# H3K9me2 upon KJ in MYC activated promoter/genes
plot <- avr1_prom[,c("dist","Me2_D","Me2_K")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")

################# Figure 7b #################
#### GCN5 correlation
data <- read.delim("siGCN5_count_enh.txt")

data$H3K9ac_log <- log2(data$siGCN5_H3K9ac.TD.Tag.Count.in.500.bp..69539510.5.Total..normalization.factor...0.14..effective.total...10000000./data$siCrtl_H3K9ac.TD.Tag.Count.in.500.bp..62677510.0.Total..normalization.factor...0.16..effective.total...10000000.)
data$BRD4_log <- log2(data$siGCN5_BRD4.TD.Tag.Count.in.500.bp..79193670.0.Total..normalization.factor...0.13..effective.total...10000000./data$siCrtl_BRD4.TD.Tag.Count.in.500.bp..61989711.0.Total..normalization.factor...0.16..effective.total...10000000.)
data <- data[data$H3K9ac_log < 10000 & data$H3K9ac_log > -100000,]
data <- data[data$BRD4_log < 10000 & data$BRD4_log > -100000,]
data <- data[complete.cases(data$BRD4_log),]
data <- data[complete.cases(data$H3K9ac_log),]


data <- data[data$siCrtl_H3K9ac.TD.Tag.Count.in.500.bp..62677510.0.Total..normalization.factor...0.16..effective.total...10000000. > 10,]

library(mltools)
std.error <- function(x) sd(x)/sqrt(length(x))


data[, "group"] <- bin_data(data$H3K9ac_log, bins=4, binType = "quantile")

bin1 <- data[data$group == "[0.00248965082931458, 1.66084000757984]",]
bin2 <- data[data$group == "[-0.344840489109252, 0.00248965082931458)",]
bin3 <- data[data$group == "[-0.78497532035298, -0.344840489109252)",]
bin4 <- data[data$group == "[-8.98299357469431, -0.78497532035298)",]

bin1 <- summarise(bin1, mean_H=mean(H3K9ac_log), sd_H=std.error(H3K9ac_log), mean_B=mean(BRD4_log), sd_B=std.error(BRD4_log))
bin2 <- summarise(bin2, mean_H=mean(H3K9ac_log), sd_H=std.error(H3K9ac_log), mean_B=mean(BRD4_log), sd_B=std.error(BRD4_log))
bin3 <- summarise(bin3, mean_H=mean(H3K9ac_log), sd_H=std.error(H3K9ac_log), mean_B=mean(BRD4_log), sd_B=std.error(BRD4_log))
bin4 <- summarise(bin4, mean_H=mean(H3K9ac_log), sd_H=std.error(H3K9ac_log), mean_B=mean(BRD4_log), sd_B=std.error(BRD4_log))

plot <- rbind(bin1,bin2,bin3,bin4)
plot$group <- c("bin1","bin2","bin3","bin4")

ggplot(data=plot, aes(y=mean_H, x = group)) +
  geom_bar(stat="identity") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_errorbar(aes(ymin=mean_H-sd_H, ymax=mean_H+sd_H, width=.1)) 

ggplot(data=plot, aes(y=mean_B, x = group)) +
  geom_bar(stat="identity") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_errorbar(aes(ymin=mean_B-sd_B, ymax=mean_B+sd_B, width=.1)) 


################# Figure 7c #################
MB3_incucyte <- read.delim("BT549_MB3_exp_R2.txt")

MB3_incucyte$DMSO <- as.numeric(gsub(",", ".", gsub("\\.", "", MB3_incucyte$DMSO)))
MB3_incucyte$MB3.50uM <- as.numeric(gsub(",", ".", gsub("\\.", "", MB3_incucyte$MB3.50uM)))
MB3_incucyte$MB3.75uM <- as.numeric(gsub(",", ".", gsub("\\.", "", MB3_incucyte$MB3.75uM)))
MB3_incucyte$MB3.100uM <- as.numeric(gsub(",", ".", gsub("\\.", "", MB3_incucyte$MB3.100uM)))
MB3_incucyte$DMSO.SD <- as.numeric(gsub(",", ".", gsub("\\.", "", MB3_incucyte$DMSO.SD)))
MB3_incucyte$MB3.50uM.SD <- as.numeric(gsub(",", ".", gsub("\\.", "", MB3_incucyte$MB3.50uM.SD)))
MB3_incucyte$MB3.75uM.SD <- as.numeric(gsub(",", ".", gsub("\\.", "", MB3_incucyte$MB3.75uM.SD)))
MB3_incucyte$MB3.100uM.SD <- as.numeric(gsub(",", ".", gsub("\\.", "", MB3_incucyte$MB3.100uM.SD)))

plot_conf <- melt(MB3_incucyte[1:5], id.vars = "Time")
plot_SD <- melt(MB3_incucyte[c("Time","DMSO.SD","MB3.50uM.SD","MB3.75uM.SD","MB3.100uM.SD")], id.vars = "Time")


plot <- cbind(plot_conf, plot_SD)
colnames(plot) <- c("Time","Name","conf","Time1","Name1","SD")

ggplot(data = plot, 
       aes(x = Time, y = conf, color = Name)) +
  geom_line(size = 1, alpha = 1) + theme(axis.title = element_text(size = 20),
                                         axis.text = element_text(size = 20, colour = "black"),
                                         panel.background = element_blank(),
                                         panel.grid = element_blank(),
                                         axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("Confluency in %") + xlab("Hours") +
  geom_errorbar(aes(ymin=conf-SD, ymax=conf+SD), width=2, size = 1) + scale_color_manual(values=c("#7F7F7F","#5A5A5A","#759D5D","#3F7A14")) + ylim(0,105)






################# Figure 7d #################
# Make bins for BT549 myc gene ranking (same as figure 2)
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
plot2 <- plot2[plot2$padj_MYC_vs_Crtl < 1,]
plot2 <- plot2[,c("symbol", "Log2FC_MYC_vs_Crtl")]
plot3 <- left_join(plot, plot2)
plot3[is.na(plot3)] = 0
plot3 <- plot3[!duplicated(plot3$symbol),]

# Only include expressed genes
plot3 <- plot3[abs(plot3$Log2FC_MYC_vs_Crtl) > 0,]

# Make GCN5 geneset
diff <- read.table("DE_analysis_MB3_BT549.txt")
diff <- diff[complete.cases(diff),]

# All significant regulated
cluster1_BT <- as.data.frame(diff[diff$padj.c1 < 0.05  ,]$symbol)
colnames(cluster1_BT) <- c("symbol")
cluster1_BT$geneset <- rep("MB3_sig")
cluster1_BT <- cluster1_BT[,c("geneset","symbol")]
clust1_tibble <- as.tibble(cluster1_BT)

# All non significant regulated
cluster2_BT <- as.data.frame(diff[diff$padj.c1 > 0.05,]$symbol)
colnames(cluster2_BT) <- c("symbol")
cluster2_BT$geneset <- rep("MB3_not")
cluster2_BT <- cluster2_BT[,c("geneset","symbol")]
clust2_tibble <- as.tibble(cluster2_BT)

# All significant down regulated
cluster3_BT <- as.data.frame(diff[diff$padj.c1 < 0.05 & diff$log2FC.c1 < 0,]$symbol)
colnames(cluster3_BT) <- c("symbol")
cluster3_BT$geneset <- rep("MB3_down")
cluster3_BT <- cluster3_BT[,c("geneset","symbol")]
clust3_tibble <- as.tibble(cluster3_BT)

# All significant up regulated
cluster4_BT <- as.data.frame(diff[diff$padj.c1 < 0.05 & diff$log2FC.c1 > 0,]$symbol)
colnames(cluster4_BT) <- c("symbol")
cluster4_BT$geneset <- rep("MB3_up")
cluster4_BT <- cluster4_BT[,c("geneset","symbol")]
clust4_tibble <- as.tibble(cluster4_BT)

# Create list of genesets for enrichment analysis
clusters <- rbind(clust1_tibble,clust2_tibble,clust3_tibble,clust4_tibble)
clusters$symbol <- as.integer(mapIds(org.Hs.eg.db, clusters$symbol,"ENTREZID", "SYMBOL"))
colnames(clusters) <- c("gs_name","entrez_gene")

# Bin gene based on MYC promoter and enhancer binding
library(mltools)
plot3[, "group"] <- bin_data(plot3$count, bins=4, binType = "quantile")

# Subset bins - bin 1 is promoter dominated and bin 4 is most enhancer domincated
bin1 <- plot3[plot3$group =="[71.57, 700.43]",]
bin2 <- plot3[plot3$group =="[25.53, 71.57)",]
bin3 <- plot3[plot3$group =="[-48.73, 25.53)",]
bin4 <- plot3[plot3$group =="[-1079.65, -48.73)",]

# Calculate enrichment for each gene bin
Genes <- (mapIds(org.Hs.eg.db, bin1$symbol,"ENTREZID", "SYMBOL"))
bin1 <- enricher(Genes, TERM2GENE=clusters, pvalueCutoff =1, maxGSSize = 10000, )
bin1 <- as.data.frame(bin1@result)
bin1$type <- rep("bin1")

Genes <- (mapIds(org.Hs.eg.db, bin2$symbol,"ENTREZID", "SYMBOL"))
bin2 <- enricher(Genes, TERM2GENE=clusters, pvalueCutoff =1, maxGSSize = 10000)
bin2 <- as.data.frame(bin2@result)
bin2$type <- rep("bin2")

Genes <- (mapIds(org.Hs.eg.db, bin3$symbol,"ENTREZID", "SYMBOL"))
bin3 <- enricher(Genes, TERM2GENE=clusters, pvalueCutoff =1, maxGSSize = 10000)
bin3 <- as.data.frame(bin3@result)
bin3$type <- rep("bin3")

Genes <- (mapIds(org.Hs.eg.db, bin4$symbol,"ENTREZID", "SYMBOL"))
bin4 <- enricher(Genes, TERM2GENE=clusters, pvalueCutoff =1, maxGSSize = 10000)
bin4 <- as.data.frame(bin4@result)
bin4$type <- rep("bin4")

# Merge for each bin
plot <- rbind(bin1,bin2,bin3,bin4)

# Subset GCN5 down-regulated genes and plot
plot_sig <- plot[plot$ID == "MB3_down",]

ggplot(data=plot_sig, aes(y=-log10(qvalue), x=ID, fill = type)) + geom_bar(stat="identity", position = "dodge") + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")  + xlab(" ") + ylab("-log10(q-value)")







agg <- read.delim("MYC_genes_all_rep.txt", header = TRUE)
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_prom)<-c("dist","V1","V2","V3","61","62","63","i1","i2","i3")

plot <- avr2_prom[,c("dist","V2","63","i3")]

plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=1)+
  scale_fill_viridis_d(option="E") + scale_color_manual(values = c("red","blue","green")) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") 





GCN5 <- read.delim("MYC_ebh_BRD4_count.txt")

GCN5$Ann <- sub(" .*","",GCN5$Annotation)
GCN5_enhancer <- GCN5[GCN5$Ann == "Intergenic" | GCN5$Ann == "intron",]

GCN5_enhancer$GCN5_LOG2_MYCMI6 <- log2(GCN5_enhancer$MYCMI6.merged.Tag.Count.in.1000.bp..66652446.0.Total..normalization.factor...0.15..effective.total...10000000./GCN5_enhancer$VEH.merged.Tag.Count.in.1000.bp..49733911.0.Total..normalization.factor...0.20..effective.total...10000000.)
GCN5_enhancer$GCN5_LOG2_MYCi975 <- log2(GCN5_enhancer$MYCi975.merged.Tag.Count.in.1000.bp..64009079.5.Total..normalization.factor...0.16..effective.total...10000000./GCN5_enhancer$VEH.merged.Tag.Count.in.1000.bp..49733911.0.Total..normalization.factor...0.20..effective.total...10000000.)
GCN5_enhancer$GCN5_LOG2_KJ <- log2(GCN5_enhancer$Combined_KJ_BRD4.TD.Tag.Count.in.1000.bp..68014560.5.Total..normalization.factor...0.15..effective.total...10000000./GCN5_enhancer$Combined_DMSO_BRD4.TD.Tag.Count.in.1000.bp..68548071.0.Total..normalization.factor...0.15..effective.total...10000000.)

GCN5_KJ <- read.delim("MYC_KJ_gene_BRD4_count.txt")

GCN5_KJ$GCN5_LOG2_KJ <- log2(GCN5_KJ$Combined_KJ_BRD4.TD.Tag.Count.in.1000.bp..68014560.5.Total..normalization.factor...0.15..effective.total...10000000./GCN5_KJ$Combined_DMSO_BRD4.TD.Tag.Count.in.1000.bp..68548071.0.Total..normalization.factor...0.15..effective.total...10000000.)
GCN5_KJ$GCN5_LOG2_MYCMI6 <- log2(GCN5_KJ$MYCMI6.merged.Tag.Count.in.1000.bp..66652446.0.Total..normalization.factor...0.15..effective.total...10000000./GCN5_KJ$VEH.merged.Tag.Count.in.1000.bp..49733911.0.Total..normalization.factor...0.20..effective.total...10000000.)
GCN5_KJ$GCN5_LOG2_MYCi975 <- log2(GCN5_KJ$MYCi975.merged.Tag.Count.in.1000.bp..64009079.5.Total..normalization.factor...0.16..effective.total...10000000./GCN5_KJ$VEH.merged.Tag.Count.in.1000.bp..49733911.0.Total..normalization.factor...0.20..effective.total...10000000.)


# Plot half violin 
ggplot() + 
  geom_boxplot(data=GCN5_enhancer, aes(y=GCN5_LOG2_KJ, x = "1_KJ_enh"),outlier.alpha = 0) +
  geom_boxplot(data=GCN5_enhancer, aes(y=GCN5_LOG2_MYCMI6, x = "2_MYCMI6_enh"),outlier.alpha = 0) +
  geom_boxplot(data=GCN5_enhancer, aes(y=GCN5_LOG2_MYCi975, x = "3_MYCi975_enh"),outlier.alpha = 0) +
  geom_boxplot(data=GCN5_KJ, aes(y=GCN5_LOG2_KJ, x = "4_KJ_prom"),outlier.alpha = 0) +
  geom_boxplot(data=GCN5_KJ, aes(y=GCN5_LOG2_MYCMI6, x = "5_MYCMI6_prom"),outlier.alpha = 0) +
  geom_boxplot(data=GCN5_KJ, aes(y=GCN5_LOG2_MYCi975, x = "MYCi975_prom"),outlier.alpha = 0) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") +
  geom_vline(xintercept = 0, linetype="dashed") + ylab("Log2Foldchange") + geom_hline(yintercept = 0, linetype="dashed") + coord_cartesian(ylim = c(-3,3))


t.test(GCN5_enhancer$GCN5_LOG2_MYCi975, GCN5_promoter$GCN5_LOG2_MYCi975, paired = F)











################# Figure 7e #################
depmap <- read.csv("CRISPR_(DepMap_Public_23Q2+Score,_Chronos)_subsetted.csv")

# Make rank plot
ggplot(data=depmap, aes(y=KAT2A, x=rank(-KAT2A), col = highlight)) + 
  geom_point(alpha = 1) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("Gene effect") + xlab("Ranked(Gene effect)") + scale_color_manual(values =c("grey","#B11518")) + geom_hline(yintercept = 0,  linetype = "dashed")


################# Supp 8a #################
annotation_T47D <- read.delim("MYC_GCN5_Common_ann.txt") # Load in dataframe
annotation_T47D$Ann <- sub(" .*","",annotation_T47D$Annotation) # Remove intron annotation
annotation_T47D_plot <- as.data.frame(table(annotation_T47D$Ann)) # Summaries data for plotting
annotation_T47D_plot$type <- rep("GCN5-MYC")

plot <- rbind(annotation_T47D_plot)

library(tidyverse)
library(viridis)

ggplot(plot, aes(x=type, y=Freq, fill  = Var1)) +  geom_bar(position="fill", stat="identity") + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "right")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("% Binding distribution") + xlab(" ") 

################# Supp 8b #################
# annotatePeaks.pl MYC_enh.bed hg38 -d tagdig/GCN5_D1.TD/ -size 5000 -hist 50 > MYC_enh_GCN5.agg.txt
Enhancer <- read.delim("MYC_enh_GCN5.agg.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","GCN5_Enh")

# annotatePeaks.pl ../CellLine_MYC_MED1_ChIP/peak/MYC_BT549_Consensus.Promoter.bed hg38 -d tagdig/GCN5_D1.TD/ -size 5000 -hist 50 > MYC_prom_GCN5.agg.txt
Promoter <- read.delim("MYC_prom_GCN5.agg.txt")
avr1_prom <-Promoter[,c(1,grep("Coverage", colnames(Promoter)))]
colnames(avr1_prom)<-c("dist","GCN5_Prom")

# Merge avarage count matrixes for promoter and enhancer bound by MYC
plot <- cbind(avr1_enh, avr1_prom)
plot <-melt(plot, id="dist")

# Make the plot
ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")
################# Supp 8e #################
# annotatePeaks.pl ../CellLine_MYC_MED1_ChIP/peak/MYC_BT549_Consensus.Promoter.bed hg38 -d tagdig/GCN5_D.TD tagdig/GCN5_K.TD -size 1000 -cpu 10 > MYC_prom_count.txt
Prom <- read.delim("MYC_prom_count.txt")

# Calculate GCN5 log2fc change upon KJ-PYR-9 treatment at ALL MYC bound promoters in BT549
Prom$log2fc <- log2(Prom$tagdig.GCN5_K.TD.Tag.Count.in.1000.bp..63138388.0.Total..normalization.factor...0.16..effective.total...10000000./Prom$tagdig.GCN5_D.TD.Tag.Count.in.1000.bp..62665990.5.Total..normalization.factor...0.16..effective.total...10000000.)

# annotatePeaks.pl MYC_enh.bed hg38 -d tagdig/GCN5_D.TD tagdig/GCN5_K.TD -size 1000 -cpu 10 > MYC_Enh_count.txt
Enh <- read.delim("MYC_Enh_count.txt")

# Calculate GCN5 log2fc change upon KJ-PYR-9 treatment at ALL MYC bound enhancers in BT549
Enh$log2fc <- log2(Enh$tagdig.GCN5_K.TD.Tag.Count.in.1000.bp..63138388.0.Total..normalization.factor...0.16..effective.total...10000000./Enh$tagdig.GCN5_D.TD.Tag.Count.in.1000.bp..62665990.5.Total..normalization.factor...0.16..effective.total...10000000.)

# Plot half violin 
ggplot() + 
  geom_flat_violin(data=Prom, aes(y=log2fc, x = "Prom")) + geom_boxplot(data=Prom, aes(y=log2fc, x = "Prom"), width=0.1, outlier.alpha = 0) +
  geom_flat_violin(data=Enh, aes(y=log2fc, x = "Enh")) + xlab("") + geom_boxplot(data=Enh, aes(y=log2fc, x = "Enh"), width=0.1,  outlier.alpha = 0) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") +
  geom_vline(xintercept = 0, linetype="dashed") + ylab("Log2Foldchange") + coord_flip() + geom_hline(yintercept = 0, linetype="dashed") + ylim(-1.5,1.5) 

Prom <- Prom[complete.cases(Prom$log2fc),]
Prom <- Prom[is.finite(Prom$log2fc),]
Enh <- Enh[complete.cases(Enh$log2fc),]
Enh <- Enh[is.finite(Enh$log2fc),]
t.test(Prom$log2fc,Enh$log2fc, paired = F, var.equal = T)$p.value

################# Supp 8f #################
#annotatePeaks.pl ../../MYC_Project/CellLine_MYC_MED1_ChIP/peak/MYC_BT549_Consensus.bed hg38 -d ../../MYC_Project/BT549_GCN_KDM_KJ/tagdig/GCN5_* ../../MYC_Project/BT549_H3K9_KJ/tagdig/H3K9ac_* *merged ../BT549_H3K9ac_inh/*merged -size 1000 > GCN5_count_myc_sites.txt
GCN5 <- read.delim("GCN5_count_myc_sites.txt")

GCN5$Ann <- sub(" .*","",GCN5$Annotation)
GCN5_enhancer <- GCN5[GCN5$Ann == "Intergenic" | GCN5$Ann == "intron",]
GCN5_promoter <- GCN5[GCN5$Ann == "promoter-TSS",]

GCN5_enhancer$GCN5_LOG2_MYCi975 <- log2(GCN5_enhancer$MYCi975.merged.Tag.Count.in.1000.bp..48925430.5.Total..normalization.factor...0.20..effective.total...10000000./GCN5_enhancer$VEH.merged.Tag.Count.in.1000.bp..57950838.0.Total..normalization.factor...0.17..effective.total...10000000.)
GCN5_promoter$GCN5_LOG2_MYCi975 <- log2(GCN5_promoter$MYCi975.merged.Tag.Count.in.1000.bp..48925430.5.Total..normalization.factor...0.20..effective.total...10000000./GCN5_promoter$VEH.merged.Tag.Count.in.1000.bp..57950838.0.Total..normalization.factor...0.17..effective.total...10000000.)
GCN5_enhancer <- GCN5_enhancer[GCN5_enhancer$GCN5_LOG2_MYCi975 < 1000 & GCN5_enhancer$GCN5_LOG2_MYCi975 > -1000,]
GCN5_promoter <- GCN5_promoter[GCN5_promoter$GCN5_LOG2_MYCi975 < 1000 & GCN5_promoter$GCN5_LOG2_MYCi975 > -1000,]

# Plot half violin 
ggplot() + 
  geom_flat_violin(data=GCN5_promoter, aes(y=GCN5_LOG2_MYCi975, x = "Prom")) + geom_boxplot(data=GCN5_promoter, aes(y=GCN5_LOG2_MYCi975, x = "Prom"), width=0.1, outlier.alpha = 0) +
  geom_flat_violin(data=GCN5_enhancer, aes(y=GCN5_LOG2_MYCi975, x = "Enh")) + xlab("") + geom_boxplot(data=GCN5_enhancer, aes(y=GCN5_LOG2_MYCi975, x = "Enh"), width=0.1,  outlier.alpha = 0) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") +
  geom_vline(xintercept = 0, linetype="dashed") + ylab("Log2Foldchange") + coord_flip() + geom_hline(yintercept = 0, linetype="dashed") + ylim(-1.5,1.5) 


t.test(GCN5_enhancer$GCN5_LOG2_MYCi975, GCN5_promoter$GCN5_LOG2_MYCi975, paired = F, var.equal = T)$p.value

################# Supp 8g #################
Enhancer <- read.delim("NewInhib_GCN5_H3K9ac_KJprom_agg.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","GCN5_MYCMI6","GCN5_MYCi975","GCN5_Veh","H3K9ac_MYCMI6","H3K9ac_MYCi975","H3K9ac_Veh")

plot <- avr1_enh[,c("dist","H3K9ac_Veh","H3K9ac_MYCi975")]
plot <-melt(plot, id="dist")

# Make the plot
ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")


plot <- avr1_enh[,c("dist","GCN5_Veh","GCN5_MYCi975")]
plot <-melt(plot, id="dist")

# Make the plot
ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")



Enhancer <- read.delim("NewInhib_GCN5_H3K9ac_agg.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","GCN5_MYCMI6","GCN5_MYCi975","GCN5_Veh","H3K9ac_MYCMI6","H3K9ac_MYCi975","H3K9ac_Veh")

plot <- avr1_enh[,c("dist","H3K9ac_Veh","H3K9ac_MYCi975")]
plot <-melt(plot, id="dist")

# Make the plot
ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")


plot <- avr1_enh[,c("dist","GCN5_Veh","GCN5_MYCi975")]
plot <-melt(plot, id="dist")

# Make the plot
ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp")






################# Supp 8h #################
# Count KDM3A binding in the four different subset of genomic regions using homer and laod into R

# annotatePeaks.pl KDM3A_only_enhancer.bed hg38 -d tagdig/KDM3A_D.TD/ tagdig/KDM3A_K.TD/ -size 1000 -cpu 12 > KDM3A_only_enhancer_count.txt
KDM3A_enh <- read.delim("KDM3A_only_enhancer_count.txt")

# annotatePeaks.pl KDM3A_only_promoter.bed hg38 -d tagdig/KDM3A_D.TD/ tagdig/KDM3A_K.TD/ -size 1000 -cpu 12 > KDM3A_only_promoter_count.txt
KDM3A_prom <- read.delim("KDM3A_only_promoter_count.txt")

# annotatePeaks.pl KDM3A_MYC_enhancer.bed hg38 -d tagdig/KDM3A_D.TD/ tagdig/KDM3A_K.TD/ -size 1000 -cpu 12 > KDM3A_MYC_share_enhancer_count.txt
MYC_KDM3A_enh <- read.delim("KDM3A_MYC_share_enhancer_count.txt")

# annotatePeaks.pl KDM3A_MYC_promoter.bed hg38 -d tagdig/KDM3A_D.TD/ tagdig/KDM3A_K.TD/ -size 1000 -cpu 12 > KDM3A_MYC_share_promoter_count.txt
MYC_KDM3A_prom <- read.delim("KDM3A_MYC_share_promoter_count.txt")

# Make the plot of log2fc of KDM3A binding at each location and represent in boxplot
ggplot() +
  geom_boxplot(data=KDM3A_enh, aes(y=log2(KDM3A_enh$tagdig.KDM3A_K.TD..Tag.Count.in.1000.bp..78212227.5.Total..normalization.factor...0.13..effective.total...10000000./KDM3A_enh$tagdig.KDM3A_D.TD..Tag.Count.in.1000.bp..59733611.5.Total..normalization.factor...0.17..effective.total...10000000.), x = "KDM3A_Enh"), outlier.alpha = 0) +
  geom_boxplot(data=MYC_KDM3A_enh, aes(y=log2(MYC_KDM3A_enh$tagdig.KDM3A_K.TD..Tag.Count.in.1000.bp..78212227.5.Total..normalization.factor...0.13..effective.total...10000000./MYC_KDM3A_enh$tagdig.KDM3A_D.TD..Tag.Count.in.1000.bp..59733611.5.Total..normalization.factor...0.17..effective.total...10000000.), x = "MYC_KDM3A_Enh"), outlier.alpha = 0) +
  geom_boxplot(data=KDM3A_prom, aes(y=log2(KDM3A_prom$tagdig.KDM3A_K.TD..Tag.Count.in.1000.bp..78212227.5.Total..normalization.factor...0.13..effective.total...10000000./KDM3A_prom$tagdig.KDM3A_D.TD..Tag.Count.in.1000.bp..59733611.5.Total..normalization.factor...0.17..effective.total...10000000.), x = "KDM3A_Prom"), outlier.alpha = 0) +
  geom_boxplot(data=MYC_KDM3A_prom, aes(y=log2(MYC_KDM3A_prom$tagdig.KDM3A_K.TD..Tag.Count.in.1000.bp..78212227.5.Total..normalization.factor...0.13..effective.total...10000000./MYC_KDM3A_prom$tagdig.KDM3A_D.TD..Tag.Count.in.1000.bp..59733611.5.Total..normalization.factor...0.17..effective.total...10000000.), x = "MYC_KDM3A_Prom"), outlier.alpha = 0) + coord_cartesian(ylim=c(-3,1)) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "right")+
  geom_hline(yintercept = 0, linetype="dashed") + ylab("KDM3A log2(KJ-PYR-9/VEH") + xlab(" ") 


KDM3A_enh <- log2(KDM3A_enh$tagdig.KDM3A_K.TD..Tag.Count.in.1000.bp..78212227.5.Total..normalization.factor...0.13..effective.total...10000000./KDM3A_enh$tagdig.KDM3A_D.TD..Tag.Count.in.1000.bp..59733611.5.Total..normalization.factor...0.17..effective.total...10000000.)
MYC_KDM3A_enh <- log2(MYC_KDM3A_enh$tagdig.KDM3A_K.TD..Tag.Count.in.1000.bp..78212227.5.Total..normalization.factor...0.13..effective.total...10000000./MYC_KDM3A_enh$tagdig.KDM3A_D.TD..Tag.Count.in.1000.bp..59733611.5.Total..normalization.factor...0.17..effective.total...10000000.)

KDM3A_enh <- KDM3A_enh[is.finite(KDM3A_enh)]
KDM3A_enh <- KDM3A_enh[complete.cases(KDM3A_enh)]
MYC_KDM3A_enh <- MYC_KDM3A_enh[is.finite(MYC_KDM3A_enh)]
MYC_KDM3A_enh <- MYC_KDM3A_enh[complete.cases(MYC_KDM3A_enh)]

t.test(KDM3A_enh, MYC_KDM3A_enh, paired = F, var.equal = T)$p.value

################# Supp 8k ################# 
siGCN5_incucyte <- read.delim("BT549_siGCN5_R2.txt")

siGCN5_incucyte$siCtl <- as.numeric(gsub(",", ".", gsub("\\.", "", siGCN5_incucyte$siCtl)))
siGCN5_incucyte$siGCN5 <- as.numeric(gsub(",", ".", gsub("\\.", "", siGCN5_incucyte$siGCN5)))
siGCN5_incucyte$ciCtl_SD <- as.numeric(gsub(",", ".", gsub("\\.", "", siGCN5_incucyte$ciCtl_SD)))
siGCN5_incucyte$siGCN5_SD <- as.numeric(gsub(",", ".", gsub("\\.", "", siGCN5_incucyte$siGCN5_SD)))

plot_conf <- melt(siGCN5_incucyte[1:3], id.vars = "Time")
plot_SD <- melt(siGCN5_incucyte[c("Time","ciCtl_SD","siGCN5_SD")], id.vars = "Time")


plot <- cbind(plot_conf, plot_SD)
colnames(plot) <- c("Time","Name","conf","Time1","Name1","SD")

ggplot(data = plot, 
       aes(x = Time, y = conf, color = Name)) +
  geom_line(size = 1, alpha = 1) + theme(axis.title = element_text(size = 20),
                                         axis.text = element_text(size = 20, colour = "black"),
                                         panel.background = element_blank(),
                                         panel.grid = element_blank(),
                                         axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("Confluency in %") + xlab("Hours") +
  geom_errorbar(aes(ymin=conf-SD, ymax=conf+SD), width=2, size = 1) + scale_color_manual(values=c("#7F7F7F","#3F7A14")) + ylim(0,105)






























