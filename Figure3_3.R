##############################################
################## FIGURE 3 ##################
##############################################

################# Figure 3A ################## 
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
a<-read.delim("MYC_BT549_Consensus.promoter.hist.txt",h=T) 
a<-read.delim("MYC_BT549_Consensus.enhancer.hist.txt",h=T) 
a2<-a[order(a$X0, decreasing = T),2:ncol(a)] 
MYC_DMSO <-a2[,1:101]
MYC_KJ <-a2[,102:202]

MYC_DMSO_hist <- Heatmap(as.matrix(MYC_DMSO), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,20,30),c("white","green1","green4")), use_raster = T, raster_quality = 10, column_split = c(rep("MYC_DMSO", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MYC_KJ_hist <- Heatmap(as.matrix(MYC_KJ), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,20,30),c("white","green1","green4")), use_raster = T, raster_quality = 10, column_split = c(rep("MYC_KJ", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))

draw(MYC_DMSO_hist+MYC_KJ_hist, heatmap_legend_side="bottom",row_title=paste(nrow(a2),"sites", collapse = " "))

agg <- read.delim("MYC_BT549_Consensus_agg_kj.txt", header = TRUE)
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


################# Figure 3B ################## 
load("test.image.Rdata")
library(DiffBind)

dba.plotPCA(bt.cons.count,  attributes=DBA_CONDITION, label=DBA_ID)


bt.cons.count <- dba.contrast(bt.cons.count, categories=DBA_CONDITION, minMembers = 2)
bt.cons.count <- dba.analyze(bt.cons.count, method=DBA_DESEQ2, bBlacklist=F, bGreylist=F)

dba.plotMA(bt.cons.count, contrast = 1, method=DBA_DESEQ2, bUsePval = 0.05)

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

write.table(diff_intergenic, "MYC_regulated_enhancers_BT549.txt",row.names = F, col.names = T, sep = "\t", quote = F)

ggplot(data=diff_intergenic, aes(x=MYC_Fold, y=-log(MYC_FDR), col=diffexpressed)) + 
  geom_point(alpha = 1)+ scale_color_manual(values = c("#B11518","grey","#B11518")) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("-log10(padj)") + xlab("log2(KJ/DMSO)") + ylim(0,50) + xlim(-10,10)

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
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("-log10(padj)") + xlab("log2(KJ/DMSO)") + ylim(0,50) + xlim(-10,10)

table(diff_promoter$diffexpressed)

write.table(diff_intergenic[diff_intergenic$diffexpressed == "UP",1:4], "MED1_MYC_Enh_Up.bed", sep = "\t", quote = F, row.names = F)
write.table(diff_promoter[diff_promoter$diffexpressed == "DOWN",2:4], "MED1_MYC_Prom_Down.bed", sep = "\t", quote = F, row.names = F)
write.table(diff_promoter[diff_promoter$diffexpressed == "NO",2:4], "MED1_MYC_Prom_No.bed", sep = "\t", quote = F, row.names = F)


################# Figure 3C ################## 
load("MCF7.image.Rdata")
library(DiffBind)

dba.plotPCA(bt.cons.count,  attributes=DBA_CONDITION, label=DBA_ID)


bt.cons.count <- dba.contrast(bt.cons.count, categories=DBA_CONDITION, minMembers = 2)
bt.cons.count <- dba.analyze(bt.cons.count, method=DBA_DESEQ2, bBlacklist=F, bGreylist=F)

dba.plotMA(bt.cons.count, contrast = 1, method=DBA_DESEQ2, bUsePval = 0.05)

MYC_report<-dba.report(bt.cons.count, th=1, DataType = DBA_DATA_FRAME, method = DBA_DESEQ2, contrast = 1)
MYC_rpkm<-bt.cons.count$peaks[[1]] #Get additional count info for each peak
MYC_rpkm$id<-rownames(MYC_rpkm)
MYC_report$id<-rownames(MYC_report)
MYC_data<-merge(MYC_report, MYC_rpkm, by="id")
colnames(MYC_data)[2:ncol(MYC_data)]<-paste("MYC", colnames(MYC_data)[2:ncol(MYC_data)],sep = "_") #Add MYC to all column names

# Add homer annotation
bed <- MYC_data[,2:4]
bed$id <- MYC_data$id
write.table(bed,"MCF7_MYC_sites.bed", quote=F, append=F, row.names = F, col.names = F, sep = "\t")
ann <- read.delim("MCF7_MYC_sites_ann.txt")
MYC_data_ann <- merge(ann, MYC_data, by.x ="PeakID..cmd.annotatePeaks.pl.MCF7_MYC_sites.bed.hg38.", by.y="id")

# Volcano plot of MCF7 #
MYC_data_ann$Ann <- sub(" .*","",MYC_data_ann$Annotation)
diff_intergenic <- MYC_data_ann[MYC_data_ann$Ann == "Intergenic" | MYC_data_ann$Ann == "intron" ,]
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
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("-log10(padj)") + xlab("log2(KJ/DMSO)") + ylim(0,30) + xlim(-10,10)

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
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("-log10(padj)") + xlab("log2(KJ/DMSO)") + ylim(0,30) + xlim(-10,10)

table(diff_promoter$diffexpressed)

write.table(diff_intergenic[diff_intergenic$diffexpressed == "UP",1:4], "MED1_MYC_Enh_Up.bed", sep = "\t", quote = F, row.names = F)
write.table(diff_intergenic[diff_intergenic$diffexpressed == "DOWN",1:4], "test_MED1_MYC_Enh_Down.bed", sep = "\t", quote = F, row.names = F)
write.table(diff_intergenic[diff_intergenic$diffexpressed == "NO",2:4], "MED1_MYC_Enh_No.bed", sep = "\t", quote = F, row.names = F)

################# Figure 3D ################## 
# mergePeaks MED1_MYC_Enh_Up.bed ../MCF7_KJ_MED1/MCF7_MED1_MYC_Enh_Up.bed
Percent <- c(5.27,18.02)
Type <- c("Enhancer","Promoter")

plot <- as.data.frame(rbind(Percent, Type))

ggplot(plot, aes(y=Percent, x=Type, fill = Type)) + geom_bar(stat="summary") + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "NA") + xlab("") + ylab("Percent overlap") + scale_fill_viridis_d(option="E")


################# Figure 3E ################# 
DEG_down <- read.delim("Cancer_MED1_MYC_Enh_Down_count.txt")
DEG_down <- DEG_down[,20:41]
colnames(DEG_down) <- gsub("\\..*","",colnames(DEG_down))
data <- as.data.frame(colMeans(DEG_down))
colnames(data) <- c("Deg")


DEG_no <- read.delim("Cancer_MED1_MYC_Enh_no_count.txt")
DEG_no <- DEG_no[,20:41]
colnames(DEG_no) <- gsub("\\..*","",colnames(DEG_no))
data$No <- (colMeans(DEG_no))


null <- c(0,0)
data <- rbind(data,null)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

data$scale_DEG <- range01(data$Deg)
data$scale_No <- range01(data$No)
data$type <- rownames(data)
data <- data[1:22,]

ggplot() + 
  geom_point(data = data, mapping=aes(y=scale_DEG, fill=type,x=reorder(type, -scale_DEG+-scale_No)), col = "#B11518", size = 4) +
  geom_point(data = data, mapping=aes(y=scale_No, fill=type,x=reorder(type, -scale_DEG+-scale_No)), col = "#979797", size = 4) +
  geom_line(data = data, mapping=aes(y=scale_DEG,x=reorder(type, -scale_DEG+-scale_No)), col = "#B11518", size = 2, group = 1) +
  geom_line(data = data, mapping=aes(y=scale_No,x=reorder(type, -scale_DEG+-scale_No)), col = "#979797", size = 2, group = 2) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "NA",  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Sccaled H3K27ac TagCounts") + xlab(" ") + ylim(0.20,1)


DEG_down <- read.delim("Cancer_MED1_MYC_Prom_Down_count.txt")
DEG_down <- DEG_down[,20:41]
colnames(DEG_down) <- gsub("\\..*","",colnames(DEG_down))
data <- as.data.frame(colMeans(DEG_down))
colnames(data) <- c("Deg")

DEG_no <- read.delim("Cancer_MED1_MYC_Prom_no_count.txt")
DEG_no <- DEG_no[,20:41]
colnames(DEG_no) <- gsub("\\..*","",colnames(DEG_no))
data$No <- (colMeans(DEG_no))


null <- c(0,0)
data <- rbind(data,null)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

data$scale_DEG <- range01(data$Deg)
data$scale_No <- range01(data$No)
data$type <- rownames(data)
data <- data[1:22,]

ggplot() + 
  geom_point(data = data, mapping=aes(y=scale_DEG, fill=type,x=reorder(type, -scale_DEG+-scale_No)), col = "#253884", size = 4) +
  geom_point(data = data, mapping=aes(y=scale_No, fill=type,x=reorder(type, -scale_DEG+-scale_No)), col = "#979797", size = 4) +
  geom_line(data = data, mapping=aes(y=scale_DEG,x=reorder(type, -scale_DEG+-scale_No)), col = "#253884", size = 2, group = 1) +
  geom_line(data = data, mapping=aes(y=scale_No,x=reorder(type, -scale_DEG+-scale_No)), col = "#979797", size = 2, group = 2) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "NA",  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Sccaled H3K27ac TagCounts") + xlab(" ") + ylim(0.20,1)








################# Figure 3F ################# 
## Giggle analysis BT549 #
Giggle_up <- read.csv("Giggle_up.csv")
Giggle_down <- read.csv("Giggle_down.csv")

ggplot(Giggle_down, aes(y=GIGGLE_score, x=reorder(Factor,-GIGGLE_score), fill = Factor)) + geom_bar(stat="summary") + geom_point() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "NA") + xlab("") + ylab("GIGGLE score") + scale_fill_viridis_d(option="E") + ylim(0,2000)

ggplot(Giggle_up, aes(y=GIGGLE_score, x=reorder(Factor,-GIGGLE_score), fill = Factor)) + geom_bar(stat="summary") + geom_point() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "NA") + xlab("") + ylab("GIGGLE score") + scale_fill_viridis_d(option="E") + ylim(0,2000)


################# Figure 3G ################# 
MYC_Enh_foot <- read.delim("bindetect_results.txt")

MYC_Enh_foot$diffexpressed <- "NO"
MYC_Enh_foot$diffexpressed[MYC_Enh_foot$DMSO_KJ_pvalue < 1e-50  & MYC_Enh_foot$DMSO_KJ_change < -0.1] <- "UP"
MYC_Enh_foot$diffexpressed[MYC_Enh_foot$DMSO_KJ_pvalue < 1e-50 & MYC_Enh_foot$DMSO_KJ_change > 0.1] <- "DOWN"

highlight<-(c("STAT3","TEAD2","FOXC1","SOX13")) 


ggplot(data=MYC_Enh_foot, aes(x=-DMSO_KJ_change, y=-log10(DMSO_KJ_pvalue), col=diffexpressed)) + 
  geom_point(alpha = 1)+ scale_color_manual(values = c("#253884","grey","#253884")) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("-log10(padj)") + xlab("Binding Score") +
        geom_text_repel(data=MYC_Enh_foot[MYC_Enh_foot$name %in% highlight,],aes(x = -DMSO_KJ_change, y = -log10(DMSO_KJ_pvalue),label=name, fill = "red"), 
                  box.padding = 2,size=10,
                  col = "black")



################# Figure 3H ################# 

protName<-"MED1" #Set genesymbol for bait protein
protID<-"Q15648" #Set uniport id for bait protein
a<-read.delim("QEHF3_094_RIME19_20210602_TR_peptides.txt",h=T) #Column names for each sample should match the metadata
a<- a[,1:15]
meta<-read.delim("metafil.txt",h=T, stringsAsFactors = F) #Import metadata
n<-9 #Number of samples in TMT experiment
fasta<-"MED1.fasta.txt" #Path to fasta file for bait protein
scaleFunc<-median #Median scaling is not too hard on the data. This is done within groups.
con<-c(MycI_vs_Crtl = "MycI - Crtl") #Describe contrast
padj<-0.2
lfc<-0
highlight<-(c("SOX13","YAP1")) #Proteins to highlight in plots


a2<-na.omit(a) #Removes missing values in the data since they have no quantitative info
a3<-a2[a2$X..Protein.Groups==1,] #Only retain those peptides that map uniquely to proteins
library(qPLEXanalyzer)
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
pcaPlot(MSnset_data, labelColumn = "SampleName", pointsize = 2,labelsize = 3) + xlim(-600,200) #PCA plot




#Scale peptide intensities according to group median (DEX/Eth vs IgG)
MSnset_norm_gs <- groupScaling(MSnset_data, 
                               groupingColumn = "Grp")

intensityPlot(MSnset_norm_gs, title = "Peptide intensity distribution")
intensityBoxplot(MSnset_norm_gs, title = "Peptide intensity distribution")
corrPlot(MSnset_norm_gs) 
hierarchicalPlot(MSnset_norm_gs) #Tree plot
pcaPlot(MSnset_norm_gs, labelColumn = "BioRep", pointsize = 5,labelsize = 0)




#Retrieve annotation for conversion of peptide info to protein info - NOT normalized data
library(UniProt.ws)
library(dplyr)
proteins <- unique(fData(MSnset_norm_gs)$Accessions)
columns <- c("ENTRY-NAME", "PROTEIN-NAMES", "GENES")
hs <- UniProt.ws::UniProt.ws(taxId = 9606)
hs_anno <- UniProt.ws::select(hs, proteins, columns, "UNIPROTKB") %>%
  as_tibble() %>%
  mutate(GeneSymbol = gsub(" .*", "", GENES)) %>%
  dplyr::select(
    Accessions = "UNIPROTKB", Gene = "ENTRY-NAME",
    Description = "PROTEIN-NAMES", GeneSymbol
  )

#Summarise peptide intensities for proteins
MSnset_Pnorm <- summarizeIntensities(MSnset_norm_gs, mean, hs_anno)

#Get diff results
diffstats <- computeDiffStats(MSnset_Pnorm, contrasts = con)

diffexp <- getContrastResults(
  diffstats = diffstats,
  contrast = "MycI - Crtl",
  controlGroup = "iIgG"
)

maVolPlot(diffstats, contrast = con, plotType = "MA", title = con)
maVolPlot(diffstats,
          contrast = con,
          plotType = "Volcano",
          title = con)


#Make final volcano plot


spec<-diffexp[diffexp$controlLogFoldChange>2,]
high<-spec[spec$GeneSymbol %in% highlight,c("Description","GeneSymbol")]
colnames(high)[2]<-"text"
spec2<-merge(spec,high, by="Description", all.x = T)
spec2$text[is.na(spec2$text)]<-""

highlight<-(c("STAT3")) 
spec2$highlight <- "NO"
spec2$highlight[spec2$log2FC > 0 & spec2$adj.P.Val < 0.1] <- "YES"
spec2$highlight[spec2$log2FC < 0 & spec2$adj.P.Val < 0.1] <- "YES"


ggplot(data=spec2, aes(x=log2FC, y=-log10(adj.P.Val), col = highlight)) + 
  geom_point(alpha = 1) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("-log10(padj)") + xlab("log2(KJ/DMSO)") +
  geom_text_repel(data=spec2[spec2$GeneSymbol %in% highlight,],aes(x = log2FC, y = -log10(adj.P.Val),label=GeneSymbol), 
                  box.padding = 1,size=10,
                  col = "black") + xlim(-1.5,1.5) + scale_color_manual(values =c("grey","#B11518"))



################# Figure 3I ################# 
STAT3_en <- read.delim("STAT3_MA0144.2_enh.txt")
STAT3_pr <- read.delim("STAT3_MA0144.2_prom.txt")

t.test(STAT3_en$DMSO_KJ_log2fc, STAT3_pr$DMSO_KJ_log2fc)

ggplot() + 
  geom_boxplot(data=STAT3_en, aes(y=-DMSO_KJ_log2fc, x = "STAT3_en"), outlier.alpha = 0) +
  geom_boxplot(data=STAT3_pr, aes(y=-DMSO_KJ_log2fc, x = "STAT3_pr"), outlier.alpha = 0) +  xlab("") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") +
  geom_vline(xintercept = 0, linetype="dashed") + ylab("Log2Foldchange") + geom_hline(yintercept = 0, linetype="dashed")+ coord_cartesian(ylim=c(-1.5,1.5))

STAT3_en$diffexpressed <- "NO"
STAT3_en$diffexpressed[STAT3_en$DMSO_KJ_log2fc > 1 ] <- "UP"
STAT3_en$diffexpressed[STAT3_en$DMSO_KJ_log2fc < -1 ] <- "DOWN"

ggplot(data=STAT3_en, aes(y=-DMSO_KJ_log2fc, x=rank(DMSO_KJ_log2fc), col=diffexpressed)) + 
  geom_point(alpha = 1)+ scale_color_manual(values = c("#253884","grey","#253884")) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("-log10(padj)") + xlab("Binding Score")
 
 
write.table(STAT3_en[STAT3_en$diffexpressed=="UP",1:3], "STAT3_dep_MYC_enh.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(STAT3_en[STAT3_en$diffexpressed=="DOWN",1:3], "STAT3_non_MYC_enh.bed", quote = F, col.names = F, row.names = F, sep = "\t")

################# Figure 3J ################# 
# TOBIAS PlotAggregate --TFBS STAT3_dep_MYC_enh.bed --signals Merged_DMSO_corrected.bw Merged_KJ_corrected.bw --output STAT3_dep_MYC.pdf --share_y both --plot_boundaries --signal-on-x --smooth 10 --flank 120

################# Figure 3K ################# 
agg <- read.delim("STAT3_dep_MYC_enh_avg.txt", header = TRUE)
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_prom)<-c("dist","MED1_DMSO","MED1_KJ")

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







#### ATAC-seq heatmaps ####
a<-read.delim("MED1_MYC_Enh_Down_ATAC_MED1_hist.txt",h=T) 
a2<-a[order(a$X0, decreasing = T),2:ncol(a)] 
ATAC_DMSO <-a2[,1:101]
ATAC_KJ <-a2[,102:202]
MED1_DMSO <-a2[,203:303]
MED1_KJ <-a2[,304:404]

ATAC_DMSO <- Heatmap(as.matrix(ATAC_DMSO), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","grey1","grey4")), use_raster = T, raster_quality = 10, column_split = c(rep("ATAC_DMSO", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(12,"cm"))
ATAC_KJ <- Heatmap(as.matrix(ATAC_KJ), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","grey1","grey4")), use_raster = T, raster_quality = 10, column_split = c(rep("ATAC_KJ", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(12,"cm"))
MED1_DMSO <- Heatmap(as.matrix(MED1_DMSO), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MED1_DMSO", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(12,"cm"))
MED1_KJ <- Heatmap(as.matrix(MED1_KJ), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MED1_KJ", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(12,"cm"))

draw(ATAC_DMSO+ATAC_KJ+MED1_DMSO+MED1_KJ, heatmap_legend_side="bottom",row_title=paste(nrow(a2),"sites", collapse = " "))


a<-read.delim("MED1_MYC_Enh_Up_ATAC_MED1_hist.txt",h=T) 
a2<-a[order(a$X0, decreasing = T),2:ncol(a)] 
ATAC_DMSO <-a2[,1:101]
ATAC_KJ <-a2[,102:202]
MED1_DMSO <-a2[,203:303]
MED1_KJ <-a2[,304:404]

ATAC_DMSO <- Heatmap(as.matrix(ATAC_DMSO), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","grey1","grey4")), use_raster = T, raster_quality = 10, column_split = c(rep("ATAC_DMSO", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(12,"cm"))
ATAC_KJ <- Heatmap(as.matrix(ATAC_KJ), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","grey1","grey4")), use_raster = T, raster_quality = 10, column_split = c(rep("ATAC_KJ", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(12,"cm"))
MED1_DMSO <- Heatmap(as.matrix(MED1_DMSO), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MED1_DMSO", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(12,"cm"))
MED1_KJ <- Heatmap(as.matrix(MED1_KJ), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,10,20),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MED1_KJ", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(12,"cm"))

draw(ATAC_DMSO+ATAC_KJ+MED1_DMSO+MED1_KJ, heatmap_legend_side="bottom",row_title=paste(nrow(a2),"sites", collapse = " "))


agg <- read.delim("MED1_MYC_Enh_Down_agg.txt", header = TRUE)
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_prom)<-c("dist","MED1_DMSO","MED1_KJ","ATAC_DMSO","ATAC_KJ")


avr2_prom <- avr2_prom[,c("dist","ATAC_DMSO","ATAC_KJ")]
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