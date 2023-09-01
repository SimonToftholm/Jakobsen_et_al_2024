##############################################
################## FIGURE 7 ##################
##############################################

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

################# Figure 7A #################
# MYC activated enhancers
agg <- read.delim("BRD4_KJ_Enh_hist.txt", header = TRUE)
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_prom)<-c("dist","DMSO","KJ")

plot <-melt(avr2_prom, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=1)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts")

# MYC activated genes
agg <- read.delim("BRD4_KJ_Prom_avg.txt", header = TRUE)
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_prom)<-c("dist","DMSO","KJ")

plot <-melt(avr2_prom, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=1)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts")




################# Figure 7B #################
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


################# Figure 7C #################
# MYC activated enhancers
ATAC <- read.delim("BT549_ATAC_KJ_JQ1_Count.txt")

ATAC$BRD_log <- log(ATAC$Combined_BRD4_JQ1.TD.Tag.Count.in.500.bp..62532442.5.Total..normalization.factor...0.16..effective.total...10000000./ATAC$Combined_BRD4_DMSO.TD.Tag.Count.in.500.bp..68484142.5.Total..normalization.factor...0.15..effective.total...10000000.)
ATAC$S5p_log <- log(ATAC$Combined_S5p_JQ1.TD.Tag.Count.in.500.bp..40487114.5.Total..normalization.factor...0.25..effective.total...10000000./ATAC$Combined_S5p_DMSO.TD.Tag.Count.in.500.bp..60870742.5.Total..normalization.factor...0.16..effective.total...10000000.)
ATAC <- ATAC[ATAC$BRD_log > -1000 & ATAC$BRD_log < 1000 &ATAC$S5p_log > -1000 & ATAC$S5p_log < 1000,]

ATAC_plot <- ATAC[,c("BRD_log","S5p_log","Combined_S5p_DMSO.TD.Tag.Count.in.500.bp..60870742.5.Total..normalization.factor...0.16..effective.total...10000000.","Combined_BRD4_DMSO.TD.Tag.Count.in.500.bp..68484142.5.Total..normalization.factor...0.15..effective.total...10000000.")]
colnames(ATAC_plot) <- c("BRD_log","S5p_log","S5p_binding","BRD4_binding")

ggplot(ATAC_plot, aes(x=(S5p_log), y=BRD_log, col = log(S5p_binding+1))) + geom_point() + scale_colour_gradient(high="#B11518", low="grey") + geom_hline(yintercept = 0,linetype ="dashed") + geom_vline(xintercept = 0,linetype ="dashed")+ ylim(-4,4) + xlim(-5,5) +
  theme(axis.title = element_text(size = 20),
       axis.text = element_text(size = 20, colour = "black"),
       panel.background = element_blank(),
       panel.grid = element_blank(),
       axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")


ggplot(ATAC_plot, aes(x=(S5p_log), y=BRD_log, col = log(BRD4_binding+1))) + geom_point() + scale_colour_gradient(high="#253884", low="grey") + geom_hline(yintercept = 0,linetype ="dashed") + geom_vline(xintercept = 0,linetype ="dashed")+ ylim(-4,4) + xlim(-5,5) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")



################# Figure 7D #################
col<- colorRampPalette(c("blue1", "blue4"))(20)
a<-read.delim("MYC_enh_hist.txt",h=T) 
a2<-a[order(a$X0, decreasing = T),2:ncol(a)] 
BRD4_DMSO <-a2[,1:101]
S5p_DMSO <-a2[,102:202]
BRD4_JQ1 <-a2[,203:303]
S5p_JQ1 <-a2[,304:404]

BRD4_DMSO <- Heatmap(as.matrix(BRD4_DMSO), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,20,40),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("BRD4_DMSO", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
BRD4_JQ1 <- Heatmap(as.matrix(BRD4_JQ1), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,20,40),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("BRD4_JQ1", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
S5p_DMSO <- Heatmap(as.matrix(S5p_DMSO), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,15,20),c("white","red1","red4")), use_raster = T, raster_quality = 10, column_split = c(rep("S5p_DMSO", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
S5p_JQ1 <- Heatmap(as.matrix(S5p_JQ1), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,15,20),c("white","red1","red4")), use_raster = T, raster_quality = 10, column_split = c(rep("S5p_JQ1", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))

draw(BRD4_DMSO+BRD4_JQ1+S5p_DMSO+S5p_JQ1, heatmap_legend_side="bottom",row_title=paste(nrow(a2),"sites", collapse = " "))


a<-read.delim("MYC_prom_hist.txt",h=T) 
a2<-a[order(a$X0, decreasing = T),2:ncol(a)] 
BRD4_DMSO <-a2[,1:101]
S5p_DMSO <-a2[,102:202]
BRD4_JQ1 <-a2[,203:303]
S5p_JQ1 <-a2[,304:404]

BRD4_DMSO <- Heatmap(as.matrix(BRD4_DMSO), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,30,40),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("BRD4_DMSO", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
BRD4_JQ1 <- Heatmap(as.matrix(BRD4_JQ1), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,30,40),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("BRD4_JQ1", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
S5p_DMSO <- Heatmap(as.matrix(S5p_DMSO), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,15,20),c("white","red1","red4")), use_raster = T, raster_quality = 10, column_split = c(rep("S5p_DMSO", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
S5p_JQ1 <- Heatmap(as.matrix(S5p_JQ1), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,15,20),c("white","red1","red4")), use_raster = T, raster_quality = 10, column_split = c(rep("S5p_JQ1", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))

draw(BRD4_DMSO+BRD4_JQ1+S5p_DMSO+S5p_JQ1, heatmap_legend_side="bottom",row_title=paste(nrow(a2),"sites", collapse = " "))

################# Figure 7E ################# 
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





################# Figure 7F ################# 
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

################# Figure 7G ################# 
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

################# Supp 1B ################# 
# MYC activated enhancers
agg <- read.delim("MYC_act_enh_BRD4_agg.txt", header = TRUE)
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_prom)<-c("dist","Veh","MYCMI6","MYCi975","DMSO","KJ")

avr2_prom <- avr2_prom[,c("dist","Veh","MYCMI6","MYCi975")]
plot <-melt(avr2_prom, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=1)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts")

# MYC activated genes
agg <- read.delim("MYC_act_genes_BRD4_agg.txt", header = TRUE)
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_prom)<-c("dist","Veh","MYCMI6","MYCi975","DMSO","KJ")

avr2_prom <- avr2_prom[,c("dist","Veh","MYCMI6","MYCi975")]
plot <-melt(avr2_prom, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=1)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts")












