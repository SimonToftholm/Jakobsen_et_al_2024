#### MYC and BRD4 ####

################# Figure 6A ################# 
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


################# Figure 6B ################# 
#### siGCN5 and BRD4, H3K9ac and RNAPII ####

#### GCN5 correlation ####
data <- read.delim("siGCN5_count_enh.txt")

data$H3K9ac_log <- log2(data$siGCN5_H3K9ac.TD.Tag.Count.in.500.bp..69539510.5.Total..normalization.factor...0.14..effective.total...10000000./data$siCrtl_H3K9ac.TD.Tag.Count.in.500.bp..62677510.0.Total..normalization.factor...0.16..effective.total...10000000.)
data$BRD4_log <- log2(data$siGCN5_BRD4.TD.Tag.Count.in.500.bp..79193670.0.Total..normalization.factor...0.13..effective.total...10000000./data$siCrtl_BRD4.TD.Tag.Count.in.500.bp..61989711.0.Total..normalization.factor...0.16..effective.total...10000000.)
data <- data[data$H3K9ac_log < 10000 & data$H3K9ac_log > -100000,]
data <- data[data$BRD4_log < 10000 & data$BRD4_log > -100000,]
data <- data[complete.cases(data$BRD4_log),]
data <- data[complete.cases(data$H3K9ac_log),]


data <- data[data$siCrtl_H3K9ac.TD.Tag.Count.in.500.bp..62677510.0.Total..normalization.factor...0.16..effective.total...10000000. > 10,]

library(mltools)
data[, "group"] <- bin_data(-data$H3K9ac_log, bins=4, binType = "quantile")


data_plot <- data[,c("group","H3K9ac_log","BRD4_log")]

data_plot <-melt(data_plot, id="group")

ggplot(data=data_plot, aes(x=group, y=(value))) +  facet_wrap(variable ~ ., scales='free_y') + 
  geom_bar(position = "dodge", stat = "summary") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")


################# Figure 6C ################# 
#### JQ1 and BRD4 and RNAPII ####
# MYC activated enhancers
ATAC <- read.delim("BT549_ATAC_KJ_JQ1_Count.txt")

ATAC$BRD_log <- log(ATAC$Combined_BRD4_JQ1.TD.Tag.Count.in.500.bp..62532442.5.Total..normalization.factor...0.16..effective.total...10000000./ATAC$Combined_BRD4_DMSO.TD.Tag.Count.in.500.bp..68484142.5.Total..normalization.factor...0.15..effective.total...10000000.)
ATAC$S5p_log <- log(ATAC$Combined_S5p_JQ1.TD.Tag.Count.in.500.bp..40487114.5.Total..normalization.factor...0.25..effective.total...10000000./ATAC$Combined_S5p_DMSO.TD.Tag.Count.in.500.bp..60870742.5.Total..normalization.factor...0.16..effective.total...10000000.)
ATAC <- ATAC[ATAC$BRD_log > -1000 & ATAC$BRD_log < 1000 &ATAC$S5p_log > -1000 & ATAC$S5p_log < 1000,]

ATAC_plot <- ATAC[,c("BRD_log","S5p_log","Combined_S5p_DMSO.TD.Tag.Count.in.500.bp..60870742.5.Total..normalization.factor...0.16..effective.total...10000000.","Combined_BRD4_DMSO.TD.Tag.Count.in.500.bp..68484142.5.Total..normalization.factor...0.15..effective.total...10000000.")]
colnames(ATAC_plot) <- c("BRD_log","S5p_log","S5p_binding","BRD4_binding")

ggplot(ATAC_plot, aes(x=(S5p_log), y=BRD_log, col = log(BRD4_binding+1))) + geom_point() + scale_colour_gradient(high="#253884", low="grey") + geom_hline(yintercept = 0,linetype ="dashed") + geom_vline(xintercept = 0,linetype ="dashed")+ ylim(-4,4) + xlim(-5,5) +
  theme(axis.title = element_text(size = 20),
       axis.text = element_text(size = 20, colour = "black"),
       panel.background = element_blank(),
       panel.grid = element_blank(),
       axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")



################# Figure 6D ################# 
library(pheatmap)
library(viridis)
library(corrplot)
library(ComplexHeatmap)
library(circlize)
library(magick)
library(reshape2)
library(ggsci)
library(VennDiagram)
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





agg <- read.delim("JQ1_enh_hist.txt", header = TRUE)
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_prom)<-c("dist","DMSO_BRD","DMSO_S5p","JQ1_BRD","JQ1_S5p")

plot <- avr2_prom[,c("dist","DMSO_S5p","JQ1_S5p")]

plot <-melt(plot, id="dist")

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
agg <- read.delim("JQ1_Prom_hist.txt", header = TRUE)
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_prom)<-c("dist","DMSO_BRD","DMSO_S5p","JQ1_BRD","JQ1_S5p")

plot <- avr2_prom[,c("dist","DMSO_BRD","JQ1_BRD")]

plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=1)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts")



















