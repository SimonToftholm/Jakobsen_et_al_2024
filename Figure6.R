##############################################
################## FIGURE 6 ##################
##############################################

# Pre-proccessed data can be downloaded from https://zenodo.org/deposit/8323614

################# Packages #################
library(tidyverse)
library(PupillometryR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(reshape)
library(clusterProfiler)
library(ggrepel)


################# Figure 6a #################
# Read in qPLEX dataframe from figure 4A
spec2 <- read.delim("MYC_DEX_qPLEX.txt")

# Select proteins to highligh
highlight<-(c("BRD4","MED1","KAT2A","STAT3","EHMT2","EP300","MAX","BRD2")) 
spec2$highlight <- "NO"

# Make rank plot
ggplot(data=spec2, aes(y=controlLogFoldChange, x=rank(controlLogFoldChange), col = highlight)) + 
  geom_point(alpha = 1) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("log2(MYC/IgG)") + xlab("Ranked") +
  geom_text_repel(data=spec2[spec2$GeneSymbol %in% highlight,],aes(x = rank(controlLogFoldChange), y = (controlLogFoldChange),label=GeneSymbol), 
                  box.padding = 1,size=10,
                  col = "black") + scale_color_manual(values =c("grey","#B11518")) + ylim(1,3)

################# Figure 6b #################
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


################# Figure 6c #################
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

################# Figure 6d #################
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

################# Supp 7a #################
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











