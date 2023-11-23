##############################################
################## FIGURE 4 ##################
##############################################

# Pre-proccessed data can be downloaded from https://zenodo.org/deposit/8323614

################# Packages################# 
library(qPLEXanalyzer)
library(tidyverse)
library(ggrepel)
library(DiffBind)
library(pheatmap)
library(viridis)
library(corrplot)
library(ComplexHeatmap)
library(circlize)
library(magick)
library(reshape2)
library(ggsci)
library(VennDiagram)

################# Figure 4a ################## 
protName<-"MYC" #Set genesymbol for bait protein
protID<-"P01106" #Set uniport id for bait protein
a<-read.delim("QEHF3_086_RIME16_frak1to9_220920_peptides.txt") #Column names for each sample should match the metadata
a<- a[,1:17]
a <- as.numeric(a)
meta<-read.delim("metafil.txt",h=T, stringsAsFactors = F) #Import metadata
n<-11 #Number of samples in TMT experiment
fasta<-"P01106_MYC_FASTA.txt" #Path to fasta file for bait protein
scaleFunc<-median #Median scaling is not too hard on the data. This is done within groups.
con<-c(DEX_vs_Eth = "DEX - Eth") #Describe contrast
padj<-0.05
lfc<-0.3
highlight<-(c("KAT2A")) #Proteins to highlight in plots


a2<-na.omit(a) #Removes missing values in the data since they have no quantitative info
a3<-a2[a2$X..Protein.Groups==1,] #Only retain those peptides that map uniquely to proteins


MSnset_data <- convertToMSnset(a3,metadata = meta,indExpData = c(7:(7+n-1)), Sequences = 2, Accessions = 6)

#Coverage plot GR
coveragePlot(MSnset_data,
             ProteinID = protID, ProteinName = protName,
             fastaFile = fasta
)

#Make intensity plots
intensityPlot(MSnset_data, title = "Peptide intensity distribution")
intensityBoxplot(MSnset_data, title = "Peptide intensity distribution")
#Correlation plots
corrPlot(MSnset_data) #Make heatmap to look at correlation of samples
hierarchicalPlot(MSnset_data) #Tree plot
pcaPlot(MSnset_data, labelColumn = "BioRep", pointsize = 5,labelsize = 0) #PCA plot




#Scale peptide intensities according to group median (DEX/Eth vs IgG)
MSnset_norm_gs <- groupScaling(MSnset_data,
                               scalingFunction = mean,
                               groupingColumn = "Grp")

intensityPlot(MSnset_norm_gs, title = "Peptide intensity distribution")
intensityBoxplot(MSnset_norm_gs, title = "Peptide intensity distribution")
corrPlot(MSnset_norm_gs) 
hierarchicalPlot(MSnset_norm_gs) #Tree plot
pcaPlot(MSnset_norm_gs, labelColumn = "BioRep", pointsize = 5,labelsize = 0)

#Retrieve annotation for conversion of peptide info to protein info - NOT normalized data
data("human_anno")

#Summarise peptide intensities for proteins
MSnset_Pnorm <- summarizeIntensities(MSnset_norm_gs, sum, human_anno)

#Get diff results
library("statmod")
diffstats <- computeDiffStats(MSnset_Pnorm, contrasts = con)
diffexp <- getContrastResults(
  diffstats = diffstats,
  contrast = con,
  controlGroup = "IgG"
)

maVolPlot(diffstats, contrast = con, plotType = "MA", title = con)
maVolPlot(diffstats,
          contrast = con,
          plotType = "Volcano",
          title = con)

#Make final volcano plot
spec<-diffexp[diffexp$controlLogFoldChange>1,]
high<-spec[spec$GeneSymbol %in% highlight,c("Description","GeneSymbol")]
colnames(high)[2]<-"text"
spec2<-merge(spec,high, by="Description", all.x = T)
spec2$text[is.na(spec2$text)]<-""
library(ggplot2)
library(ggrepel)
gg1<-ggplot(spec2,aes(x=log2FC,y=-log10(P.Value)))
p<-geom_point(colour=ifelse(gg1$data$adj.P.Val<=padj & abs(gg1$data$log2FC)>=lfc, "#b20000", "#00000065"),
              size=ifelse(gg1$data$adj.P.Val<=padj & abs(gg1$data$log2FC)>=lfc, 2, 1))
axes<-theme(aspect.ratio=1,axis.line=element_line(size=0.5),
            axis.title.x=element_text(face="bold", size=14),
            axis.title.y=element_text(face="bold", size=14),
            axis.text = element_text(size=12, colour = "black"),
            plot.title=element_text(face="bold", size=16)
)
labels<-labs(
  x="Log2 fold change", y="-log10(P-value)", title=paste(paste(paste("Volcano plot - padj <", padj),", Log2FC>"),lfc)
)
scales<-scale_x_continuous(limits=c(-1.5,1.5))
l<-geom_vline(xintercept=0,color="black", size=0.5, linetype="dashed")
repel_text<-geom_text_repel(aes(label=gg1$data$text),box.padding = unit(0.25,"lines"),point.padding = unit(0.35,"lines"),force=3,size=6)
gg1+p+theme_bw()+axes+labels+scales+l+repel_text




highlight<-(c("BRD4","MED1","KAT2A","STAT3","EHMT2","EP300","MAX")) 
spec2$highlight <- "NO"


################# Figure 4b ################## 
# Diffbind analysis was performed as described in figure 3C-D
load("GR.02022022.image.Rdata")

dba.plotPCA(bt.cons.count,  attributes=DBA_CONDITION, label=DBA_ID)

bt.cons.count <- dba.contrast(bt.cons.count, categories=DBA_CONDITION, minMembers = 2)
bt.cons.count <- dba.analyze(bt.cons.count, method=DBA_DESEQ2, bBlacklist=F, bGreylist=F)


dba.plotMA(bt.cons.count, contrast = 1, method=DBA_DESEQ2)
dba.plotMA(bt.cons.count, contrast = 10, method=DBA_DESEQ2)
dba.plotMA(bt.cons.count, contrast = 15, method=DBA_DESEQ2)


GR_report<-dba.report(bt.cons.count, th=1, DataType = DBA_DATA_FRAME, method = DBA_DESEQ2, contrast = 1)
GR_rpkm<-bt.cons.count$peaks[[1]] #Get additional count info for each peak
GR_rpkm$id<-rownames(GR_rpkm)
GR_report$id<-rownames(GR_report)
GR_data<-merge(GR_report, GR_rpkm, by="id")
colnames(GR_data)[2:ncol(GR_data)]<-paste("GR", colnames(GR_data)[2:ncol(GR_data)],sep = "_") #Add GR to all column names

MYC_report<-dba.report(bt.cons.count, th=1, DataType = DBA_DATA_FRAME, method = DBA_DESEQ2, contrast = 15)
MYC_rpkm<-bt.cons.count$peaks[[7]] #Get additional count info for each peak
MYC_rpkm$id<-rownames(MYC_rpkm)
MYC_report$id<-rownames(MYC_report)
MYC_data<-merge(MYC_report, MYC_rpkm, by="id")
colnames(MYC_data)[2:ncol(MYC_data)]<-paste("MYC", colnames(MYC_data)[2:ncol(MYC_data)],sep = "_") #Add MYC to all column names

combi<-merge(GR_data, MYC_data, by="id")

# Subset sites both bound by GR and MYC and sites occupied alone of either of these
shared<-combi[combi$`GR_p-value`<0.1 & combi$GR_Fold < 0 & combi$GR_Conc_GR_DEX > 5 & combi$MYC_Conc_MYC_DEX > 5,] # Sites were MYC is already precent but GR is recruited
gained_MYC_GR<-shared[shared$`MYC_p-value`<0.1 & shared$MYC_Fold > 0  & shared$GR_Fold < 0,] # GR gained sites where MYC is also recruited
uncahnged_MYC_GR<-shared[shared$`MYC_p-value`>0.1 & shared$GR_Fold < 0,] # GR gained sites where MYC is also recruited

uncahnged_MYC_only<-combi[combi$GR_Conc_GR_DEX < 5 & combi$MYC_Conc_MYC_DEX > 5,] # GR gained sites where MYC is also recruited
uncahnged_GR_only<-combi[combi$GR_Conc_GR_DEX > 5 & combi$MYC_Conc_MYC_DEX < 5,] # GR gained sites where MYC is also recruited

# Save as bed files and annotate with homer to only get enhancers
write.table(gained_MYC_GR[,2:4],"Gained_MYC_Diff_DiffBindsites.bed", quote=F, append=F, row.names = F, col.names = F, sep = "\t")
write.table(uncahnged_MYC_GR[,2:4],"Constand_MYC_Diff_DiffBindsites.bed", quote=F, append=F, row.names = F, col.names = F, sep = "\t")
write.table(uncahnged_MYC_only[,2:4],"uncahnged_MYC_only.bed", quote=F, append=F, row.names = F, col.names = F, sep = "\t")
write.table(uncahnged_GR_only[,2:4],"uncahnged_GR_only.bed", quote=F, append=F, row.names = F, col.names = F, sep = "\t")

################# Figure 4c ################## 
# Use above bed files on cluster to generate heatmaps
#Make heatmap of the two above groups
a<-read.delim("Gained_MYC_Diff_DiffBindsites_Enhancer_GR_ghist.txt",h=T) 
a2<-a[order(a$X0.2, decreasing = T),2:ncol(a)] 
GR_DEX <-a2[,1:101]
GR_VEH <-a2[,102:202]

GR_VEH_hist <- Heatmap(as.matrix(GR_VEH), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,50,100),c("white","red1","red4")), use_raster = T, raster_quality = 10, column_split = c(rep("GR_VEH", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(9,"cm"), width = unit(1.5,"cm"))
GR_DEX_hist <- Heatmap(as.matrix(GR_DEX), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,50,100),c("white","red1","red4")), use_raster = T, raster_quality = 10, column_split = c(rep("GR_DEX", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(9,"cm"), width = unit(1.5,"cm"))

a<-read.delim("Gained_MYC_Diff_DiffBindsites_Enhancer_MYC_MED1_ghist.txt",h=T) 
a2<-a[order(a$X0.2, decreasing = T),2:ncol(a)] 
MYC_DEX <-a2[,203:303]
MYC_VEH <-a2[,304:404]
MED1_DEX <-a2[,405:505]
MED1_VEH <-a2[,506:606]

MED1_VEH_hist <- Heatmap(as.matrix(MED1_VEH), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,15,20),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MED1_VEH", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(9,"cm"), width = unit(1.5,"cm"))
MED1_DEX_hist <- Heatmap(as.matrix(MED1_DEX), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,15,20),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MED1_DEX", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(9,"cm"), width = unit(1.5,"cm"))
MYC_VEH_hist <- Heatmap(as.matrix(MYC_VEH), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","green3","green4")), use_raster = T, raster_quality = 10, column_split = c(rep("MYC_VEH", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(9,"cm"), width = unit(1.5,"cm"))
MYC_DEX_hist <- Heatmap(as.matrix(MYC_DEX), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","green3","green4")), use_raster = T, raster_quality = 10, column_split = c(rep("MYC_DEX", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(9,"cm"), width = unit(1.5,"cm"))

draw(GR_VEH_hist+GR_DEX_hist+MYC_VEH_hist+MYC_DEX_hist+MED1_VEH_hist+MED1_DEX_hist, heatmap_legend_side="bottom",row_title=paste(nrow(a2),"sites", collapse = " "))


#Average plot
agg <- read.delim("Gained_MYC_Diff_DiffBindsites_Enhancer_agg.txt", header = TRUE)
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_prom)<-c("dist","MYC_DSMO","MYC_KJ","MED1_DSMO","MED1_KJ")

library(reshape)
plot <- avr2_prom[,c("dist","MED1_DSMO","MED1_KJ")]
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


agg <- read.delim("GR_Gained_MYC_Diff_DiffBindsites_ann.agg.txt", header = TRUE)
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_prom)<-c("dist","GR_DEX","GR_VEH")

library(reshape)
plot <- avr2_prom[,c("dist","GR_DEX","GR_VEH")]
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






################# Figure 4d ################## 
col<- colorRampPalette(c("blue1", "blue4"))(20)

a<-read.delim("Gained_MYC_Diff_DiffBindsites_Enhancer_KJ_MED1_hist.txt",h=T) 
a2<-a[order(a$X0.2, decreasing = T),2:ncol(a)] 
MED1_DEX_KJ <-a2[,1:101]
MED1_DEX <-a2[,102:202]
MYC_KJ <-a2[,203:303]
MYC_VEH <-a2[,304:404]

MED1_DEX_KJ_hist <- Heatmap(as.matrix(MED1_DEX_KJ), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,20,25),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MED1_DEX_KJ", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MED1_DEX_hist <- Heatmap(as.matrix(MED1_DEX), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,20,25),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MED1_DEX", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MYC_KJ_hist <- Heatmap(as.matrix(MYC_KJ), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,20,25),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MYC_KJ", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MYC_VEH_hist <- Heatmap(as.matrix(MYC_VEH), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,20,25),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MYC_VEH", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))

draw(MYC_VEH_hist+MYC_KJ_hist+MED1_DEX_hist+MED1_DEX_KJ_hist, heatmap_legend_side="bottom",row_title=paste(nrow(a2),"sites", collapse = " "))


################# Figure 4e ################## 

# Promoter enhancer and gene expression
# Get MYC binding at promoter and enhancer in BT549
promoter <- read.delim("tss_hg38_MYC_DEX_count.txt")
enhancer <- read.delim("GR_MYC_sites_DESeq_Intergenic.count.txt")

merged <- merge(enhancer, promoter, by = "Gene.Name")

# Add gene expression ifmorantion
merged <- merged[!duplicated(merged$Start.x),]
KJ_Gene <- read.table("DE_analysis_KJ_BT549_DEX_new.txt")
merged <- merge(merged, KJ_Gene, by.x = "Gene.Name", by.y = "symbol")

# Subset promoters with high and low MYC binding
MYC_prom_low <- merged[merged$...tagdig.MYC_DEX_merged..Tag.Count.in.1000.bp..49362524.0.Total..normalization.factor...0.20..effective.total...10000000..y < 15, ]
MYC_prom_high <- merged[merged$...tagdig.MYC_DEX_merged..Tag.Count.in.1000.bp..49362524.0.Total..normalization.factor...0.20..effective.total...10000000..y > 15, ]

# Make bed files for the enhancer and promoter beregions
write.table(MYC_prom_low[,3:5],"MYC_prom_low_Enh.bed", quote=F, append=F, row.names = F, col.names = F, sep = "\t")
write.table(MYC_prom_high[,3:5],"MYC_prom_high_Enh.bed", quote=F, append=F, row.names = F, col.names = F, sep = "\t")
write.table(MYC_prom_low[,22:24],"MYC_prom_low_Prom.bed", quote=F, append=F, row.names = F, col.names = F, sep = "\t")
write.table(MYC_prom_high[,22:24],"MYC_prom_high_Prom.bed", quote=F, append=F, row.names = F, col.names = F, sep = "\t")

# Aggregate plot of the regions high/low
agg <- read.delim("MYC_prom_high_Enh.agg.txt", header = TRUE)
avr2_no <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_no)<-c("dist","MYC_high")

agg <- read.delim("MYC_prom_low_Enh.agg.txt", header = TRUE)
avr2_yes <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_yes)<-c("dist","MYC_low")

plot <- cbind( avr2_yes,avr2_no)
plot <- plot[,c("dist","MYC_low","MYC_high")]
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

agg <- read.delim("MYC_prom_high_Prom.agg.txt", header = TRUE)
avr2_no <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_no)<-c("dist","MYC_high")

agg <- read.delim("MYC_prom_low_Prom.agg.txt", header = TRUE)
avr2_yes <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_yes)<-c("dist","MYC_low")

plot <- cbind( avr2_yes,avr2_no)
plot <- plot[,c("dist","MYC_low","MYC_high")]
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

################# Figure 4f ################## 
# Genes expression
MYC_prom_low$type <- rep("low")
MYC_prom_high$type <- rep("high")

plot <- rbind(MYC_prom_low, MYC_prom_high)
plot <- plot[,c("type","Crtl_DEX1")]
plot <- melt(plot, id = "type")

# Subset genes based on MYC binding at enhancer/promoters and distance to TSS
merged_gene_sig <- merged[merged$padj_Crtl.DEX_vs_Crtl < 1,]
merged_gene_sig <- merged_gene_sig[!duplicated(merged_gene_sig$Gene.Name),]
merged_gene_sig <- merged_gene_sig[complete.cases(merged_gene_sig$Log2FC_Crtl.DEX_vs_Crtl),]
MYC_prom_low_sig <- merged_gene_sig[merged_gene_sig$...tagdig.MYC_DEX_merged..Tag.Count.in.1000.bp..49362524.0.Total..normalization.factor...0.20..effective.total...10000000..y < 15 & merged_gene_sig$Log2FC_Crtl.DEX_vs_Crtl > 0 & abs(merged_gene_sig$Distance.to.TSS.x) < 50000, ]
MYC_prom_high_sig <- merged_gene_sig[merged_gene_sig$...tagdig.MYC_DEX_merged..Tag.Count.in.1000.bp..49362524.0.Total..normalization.factor...0.20..effective.total...10000000..y > 15 & merged_gene_sig$Log2FC_Crtl.DEX_vs_Crtl > 0 & abs(merged_gene_sig$Distance.to.TSS.x) < 50000,]

# Prepare dataframe for plot
MYC_prom_low_sig$type <- rep("low")
MYC_prom_high_sig$type <- rep("high")

plot <- rbind(MYC_prom_low_sig, MYC_prom_high_sig)
plot <- plot[,c("type","Log2FC_Crtl.DEX_vs_Crtl","Log2FC_MYC_DEX_vs_MYC")]
plot <- melt(plot, id = "type")

# Make the boxplot for high and low bound MYC sites
ggplot(plot, aes(x=type, y=value, fill = variable))  + geom_boxplot(alpha = 1, outlier.alpha = 0) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  xlab(" ") + ylab("Mean Expression log Intensity") + scale_fill_viridis_d(option="E") + coord_cartesian(ylim=c(-2,3)) + geom_hline(yintercept=0, linetype="dashed", color = "black")

# Statistical testing
t.test(MYC_prom_high_sig$Log2FC_Crtl.DEX_vs_Crtl,
       MYC_prom_high_sig$Log2FC_MYC_DEX_vs_MYC, paired = T)$p.value






################# Supp 5a ################# 
a<-read.delim("Shared_DiffBindsites_GR.hist",h=T) 
a2<-a[order(a$X0, decreasing = T),2:ncol(a)] 
GR_DEX <-a2[,1:101]
GR_VEH <-a2[,102:202]

GR_VEH_hist <- Heatmap(as.matrix(GR_VEH), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,75,100),c("white","red1","red4")), use_raster = T, raster_quality = 10, column_split = c(rep("GR_VEH", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(9,"cm"), width = unit(1.5,"cm"))
GR_DEX_hist <- Heatmap(as.matrix(GR_DEX), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,75,100),c("white","red1","red4")), use_raster = T, raster_quality = 10, column_split = c(rep("GR_DEX", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(9,"cm"), width = unit(1.5,"cm"))

draw(GR_VEH_hist+GR_DEX_hist, heatmap_legend_side="bottom",row_title=paste(nrow(a2),"sites", collapse = " "))

a<-read.delim("Shared_DiffBindsites_MYC_MED1.hist",h=T) 
a2<-a[order(a$X0, decreasing = T),2:ncol(a)] 
MYC_DEX <-a2[,203:303]
MYC_VEH <-a2[,304:404]
MED1_DEX <-a2[,405:505]
MED1_VEH <-a2[,506:606]

MED1_VEH_hist <- Heatmap(as.matrix(MED1_VEH), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,30,40),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MED1_VEH", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(9,"cm"), width = unit(1.5,"cm"))
MED1_DEX_hist <- Heatmap(as.matrix(MED1_DEX), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,30,40),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MED1_DEX", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(9,"cm"), width = unit(1.5,"cm"))
MYC_VEH_hist <- Heatmap(as.matrix(MYC_VEH), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,30,40),c("white","green3","green4")), use_raster = T, raster_quality = 10, column_split = c(rep("MYC_VEH", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(9,"cm"), width = unit(1.5,"cm"))
MYC_DEX_hist <- Heatmap(as.matrix(MYC_DEX), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,30,40),c("white","green3","green4")), use_raster = T, raster_quality = 10, column_split = c(rep("MYC_DEX", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(9,"cm"), width = unit(1.5,"cm"))

draw(GR_VEH_hist+GR_DEX_hist+MED1_VEH_hist+MED1_DEX_hist+MYC_VEH_hist+MYC_DEX_hist, heatmap_legend_side="bottom",row_title=paste(nrow(a2),"sites", collapse = " "))


