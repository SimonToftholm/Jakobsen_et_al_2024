##############################################
################## FIGURE 6 ##################
##############################################

################# Figure 6A ################# 
KDM3A_enh <- read.delim("KDM3A_only_enhancer_count.txt")
KDM3A_prom <- read.delim("KDM3A_only_promoter_count.txt")
MYC_KDM3A_enh <- read.delim("KDM3A_MYC_share_enhancer_count.txt")
MYC_KDM3A_prom <- read.delim("KDM3A_MYC_share_promoter_count.txt")


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
  

################# Figure 6B #################
Enhancer <- read.delim("KDM3A_MYC_enh_Agg.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","Ac_D","Ac_HD","Ac_HK","Ac_K",
                      "Me1_D","Me1_HD","Me1_HK","Me1_K"
                      ,"Me2_D","Me2_HD","Me2_HK","Me2_K")

plot <- avr1_enh[,c("dist","Ac_D","Ac_HD","Ac_HK")]
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



Enhancer <- read.delim("KDM3A_MYC_prom_Agg.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","Me1_D","Me1_HD","Me1_HK","Me1_K"
                      ,"Me2_D","Me2_HD","Me2_HK","Me2_K")

plot <- avr1_enh[,c("dist","Me2_D","Me2_HD","Me2_HK")]
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




################# Figure 6D ################# 
Enhancer <- read.delim("MYC_enh_GCN5.agg.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","GCN5_Enh")

Promoter <- read.delim("MYC_prom_GCN5.agg.txt")
avr1_prom <-Promoter[,c(1,grep("Coverage", colnames(Promoter)))]
colnames(avr1_prom)<-c("dist","GCN5_Prom")

plot <- cbind(avr1_enh, avr1_prom)
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

################# Figure 6E ################# 
Prom <- read.delim("MYC_prom_count.txt")
Prom$log2fc <- log2(Prom$tagdig.GCN5_K.TD.Tag.Count.in.1000.bp..63138388.0.Total..normalization.factor...0.16..effective.total...10000000./Prom$tagdig.GCN5_D.TD.Tag.Count.in.1000.bp..62665990.5.Total..normalization.factor...0.16..effective.total...10000000.)
Enh <- read.delim("MYC_Enh_count.txt")
Enh$log2fc <- log2(Enh$tagdig.GCN5_K.TD.Tag.Count.in.1000.bp..63138388.0.Total..normalization.factor...0.16..effective.total...10000000./Enh$tagdig.GCN5_D.TD.Tag.Count.in.1000.bp..62665990.5.Total..normalization.factor...0.16..effective.total...10000000.)


ggplot() + 
  geom_flat_violin(data=Prom, aes(y=log2fc, x = "Prom")) + geom_boxplot(data=Prom, aes(y=log2fc, x = "Prom"), width=0.1, outlier.alpha = 0) +
  geom_flat_violin(data=Enh, aes(y=log2fc, x = "Enh")) + xlab("") + geom_boxplot(data=Enh, aes(y=log2fc, x = "Enh"), width=0.1,  outlier.alpha = 0) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") +
  geom_vline(xintercept = 0, linetype="dashed") + ylab("Log2Foldchange") + coord_flip() + geom_hline(yintercept = 0, linetype="dashed") + ylim(-1.5,1.5) 

################# Figure 6F ################# 
Enhancer <- read.delim("MED1_MYC_Enh_Down_GCN5_KDM3A_H3K9mod.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","GCN5_D","GCN5_K","KDM3A_D","KDM3A_HD","KDM3A_HK","KDM3A_K"
                      ,"Ac_D","Ac_HD","Ac_HK","Ac_K"
                      ,"Me1_D","Me1_HD","Me1_HK","Me1_K"
                      ,"Me2_D","Me2_HD","Me2_HK","Me2_K"
                      ,"Input_D","Input_HD","Input_HK","Input_K")

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


Enhancer <- read.delim("MED1_MYC_Enh_Up_GCN5_KDM3A_H3K9mod.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","GCN5_D","GCN5_K","KDM3A_D","KDM3A_HD","KDM3A_HK","KDM3A_K"
                      ,"Ac_D","Ac_HD","Ac_HK","Ac_K"
                      ,"Me1_D","Me1_HD","Me1_HK","Me1_K"
                      ,"Me2_D","Me2_HD","Me2_HK","Me2_K"
                      ,"Input_D","Input_HD","Input_HK","Input_K")


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



Enhancer <- read.delim("MED1_MYC_Gene_Down_GCN5_KDM3A_H3K9mod.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","GCN5_D","GCN5_K","KDM3A_D","KDM3A_HD","KDM3A_HK","KDM3A_K"
                      ,"Ac_D","Ac_HD","Ac_HK","Ac_K"
                      ,"Me1_D","Me1_HD","Me1_HK","Me1_K"
                      ,"Me2_D","Me2_HD","Me2_HK","Me2_K"
                      ,"Input_D","Input_HD","Input_HK","Input_K")


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



Enhancer <- read.delim("MYC_Gene_down_MED1_agg.txt")
avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","MED1_D","MED1_K")


plot <- avr1_enh[,c("dist","MED1_D","MED1_K")]
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







#### DIFFBIND ####
load("H3K9ac.image.Rdata")
load("GCN5_GCN5sites.image.Rdata")
library(DiffBind)

dba.plotPCA(bt.cons.count,  attributes=DBA_CONDITION, label=DBA_ID)


bt.cons.count <- dba.contrast(bt.cons.count, categories=DBA_CONDITION, minMembers = 2)
bt.cons.count <- dba.analyze(bt.cons.count, method=DBA_DESEQ2, bBlacklist=F, bGreylist=F)

dba.plotMA(bt.cons.count, contrast = 1, method=DBA_DESEQ2, bUsePval = 0.05)
dba.plotHeatmap(bt.cons.count)

MYC_report<-dba.report(bt.cons.count, th=1, DataType = DBA_DATA_FRAME, method = DBA_DESEQ2, contrast = 1)
MYC_rpkm<-bt.cons.count$peaks[[1]] #Get additional count info for each peak
MYC_rpkm$id<-rownames(MYC_rpkm)
MYC_report$id<-rownames(MYC_report)
MYC_data<-merge(MYC_report, MYC_rpkm, by="id")
colnames(MYC_data)[2:ncol(MYC_data)]<-paste("MYC", colnames(MYC_data)[2:ncol(MYC_data)],sep = "_") #Add MYC to all column names

# Add homer annotation
bed <- MYC_data[,2:4]
bed$id <- MYC_data$id
write.table(bed,"MYC_sites.bed", quote=F, append=F, row.names = F, col.names = F, sep = "\t")
ann <- read.delim("MYC_sites.ann.txt")
MYC_data_ann <- merge(ann, MYC_data, by.x ="PeakID..cmd.annotatePeaks.pl.MYC_sites.bed.hg38.", by.y="id")

#### Volcano plot for enhancer and promoter regions BT549
MYC_data_ann$Ann <- sub(" .*","",MYC_data_ann$Annotation)
diff_intergenic <- MYC_data_ann[MYC_data_ann$Ann == "Intergenic" | MYC_data_ann$Ann == "intron",]
diff_promoter <- MYC_data_ann[MYC_data_ann$Ann == "promoter-TSS",]

# Sig enhancer 
diff_intergenic$diffexpressed <- "NO"
diff_intergenic$diffexpressed[diff_intergenic$MYC_Fold > 0 & diff_intergenic$MYC_FDR < 0.1] <- "UP"
diff_intergenic$diffexpressed[diff_intergenic$MYC_Fold < 0 & diff_intergenic$MYC_FDR < 0.1] <- "DOWN"

ggplot(data=diff_intergenic, aes(x=MYC_Fold, y=-log(MYC_FDR), col=diffexpressed)) + 
  geom_point(alpha = 1)+ scale_color_manual(values = c("#B11518","grey","#B11518")) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("-log10(padj)") + xlab("log2(KJ/DMSO)") + ylim(0,30) + xlim(-5,5)

table(diff_intergenic$diffexpressed)

# Sig promoter 
diff_promoter$diffexpressed <- "NO"
diff_promoter$diffexpressed[diff_promoter$MYC_FDR < 0.1  & diff_promoter$MYC_Fold > 0] <- "UP"
diff_promoter$diffexpressed[diff_promoter$MYC_FDR < 0.1 & diff_promoter$MYC_Fold < 0] <- "DOWN"

ggplot(data=diff_promoter, aes(x=MYC_Fold, y=-log(MYC_FDR), col=diffexpressed)) + 
  geom_point(alpha = 1)+ scale_color_manual(values = c("#253884","grey","#253884")) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("-log10(padj)") + xlab("log2(KJ/DMSO)") + ylim(0,30) + xlim(-5,5)

table(diff_promoter$diffexpressed)

write.table(diff_intergenic[diff_intergenic$diffexpressed == "UP",1:4], "MED1_MYC_Enh_Up.bed", sep = "\t", quote = F, row.names = F)
write.table(diff_promoter[diff_promoter$diffexpressed == "DOWN",2:4], "MED1_MYC_Prom_Down.bed", sep = "\t", quote = F, row.names = F)
write.table(diff_promoter[diff_promoter$diffexpressed == "NO",2:4], "MED1_MYC_Prom_No.bed", sep = "\t", quote = F, row.names = F)


sig <- MYC_data_ann[MYC_data_ann$`MYC_p-value` < 0.05,]
