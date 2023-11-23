##############################################
################## FIGURE 3 ##################
##############################################

# Pre-proccessed data can be downloaded from https://zenodo.org/deposit/8323614

################# Packages ################# 
library(pheatmap)
library(viridis)
library(corrplot)
library(ComplexHeatmap)
library(circlize)
library(magick)
library(reshape2)
library(ggsci)
library(VennDiagram)
library(DiffBind)
library(ggrepel)
library(qPLEXanalyzer)
library(UniProt.ws)
library(dplyr)



################# Figure 3a ################## 
# Create heatmap for MYC bound enhancer and promoters in BT549 cells upon KJ-PYR-9 treatment for 3hr. 
# Run the following on cluster for enhancers and promoter MYC bed file
## annotePeaks *bed* hg38 -d BT549_DMSO_merged.TD BT549_KJ_merged.TD -raw -size 10000 -hist 100 -ghist 

# Load in dataframe for promoter or enhancer regions
a<-read.delim("MYC_BT549_Consensus.promoter.hist.txt",h=T) 
a<-read.delim("MYC_BT549_Consensus.enhancer.hist.txt",h=T) 

# Sort dataframe by decreasing signal
a2<-a[order(a$X0, decreasing = T),2:ncol(a)] 

# Subset columns for treatment to each heatmap
MYC_DMSO <-a2[,1:101]
MYC_KJ <-a2[,102:202]

# Plot the heatmaps
MYC_DMSO_hist <- Heatmap(as.matrix(MYC_DMSO), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,20,30),c("white","green1","green4")), use_raster = T, raster_quality = 10, column_split = c(rep("MYC_DMSO", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MYC_KJ_hist <- Heatmap(as.matrix(MYC_KJ), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,20,30),c("white","green1","green4")), use_raster = T, raster_quality = 10, column_split = c(rep("MYC_KJ", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))

# Draw the two heatmaps side by side
draw(MYC_DMSO_hist+MYC_KJ_hist, heatmap_legend_side="bottom",row_title=paste(nrow(a2),"sites", collapse = " "))

# Count avarage count using homer for avarage plot of all MYC regions in BT549 cells
## annotatePeaks BT549_MYC_Concensus.bed hg38 -d BT549_DMSO_merged.TD BT549_KJ_merged.TD -size 5000 -hist 50 -raw 12 -cpu 12
agg <- read.delim("MYC_BT549_Consensus_agg_kj.txt", header = TRUE)

# Grep covarage from regions
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]

# Insert correct names
colnames(avr2_prom)<-c("dist","DMSO","KJ")

# Prepare dataframe for plotting
plot <-melt(avr2_prom, id="dist")

# Plot the average plot
ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=1) +
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts")


################# Figure 3b ################## 
# Diffbind was used to identify differential used MYC enhancers by MED1 recruitment upon KJ-PYR-9 treatment
## The diffbind analysis was run on cluster as it takes quite some rescources

#library(DiffBind)

# Read in the samplesheet
#samples<-read.csv("SampleSheetDiffBind.csv", sep = ";")
#samples<- samples[1:18,1:10]
#bt<-dba(sampleSheet = samples)

# Find concensuspeaks
#bt.cons<- dba.peakset(bt, consensus=DBA_CONDITION, minOverlap=0.99)
#bt.cons <- dba(bt.cons, mask=bt.cons$masks$Consensus, minOverlap=1)
#consensus_peaks <- dba.peakset(bt.cons, bRetrieve=TRUE)

# Count reads in peaks
#bt.cons.count<-dba.count(bt, peaks = consensus_peaks)

# Save image to run the downstream analysis local
#save.image("MED1_KJ_MYC_regions.image.Rdata)

# Load in image from analysis
load("MED1_KJ_MYC_regions.image.Rdata")

# PCA plot to validate replicates
dba.plotPCA(bt.cons.count,  attributes=DBA_CONDITION, label=DBA_ID)

# Set contrast and run DEseq2 analysis
bt.cons.count <- dba.contrast(bt.cons.count, categories=DBA_CONDITION, minMembers = 2)
bt.cons.count <- dba.analyze(bt.cons.count, method=DBA_DESEQ2, bBlacklist=F, bGreylist=F)

dba.plotMA(bt.cons.count, contrast = 1, method=DBA_DESEQ2, bUsePval = 0.05)

# Create datafrom containing the information from the DESeq2 analysis
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


################# Figure 3c ################## 
# Same analysis ad in Figure 3B but for MCF7 cells with KJ-PYR-9 treatment and MED1 ChIP-seq
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

################# Figure 3d ################# 
# Count H3K27ac in MYC activated enhancers 
## annotatePeaks.pl MED1_MYC_Enh_Down.bed hg38 -d 76F2V.bam.TD A493.bam.TD AU565.bam.TD BT549.bam.TD H128.bam.TD H2171.bam.TD HCC1937.bam.TD HCC1954.bam.TD HCT116.bam.TD HeLa.bam.TD MB231.bam.TD MB361.bam.TD MB436.bam.TD MB468.bam.TD MCF7.bam.TD MM1S.bam.TD SKBR3.bam.TD T47D.bam.TD U2OS.bam.TD U87.bam.TD UACC812.bam.TD ZR751.bam.TD -size 1000 -cpu 24

# Load into R
DEG_down <- read.delim("Cancer_MED1_MYC_Enh_Down_count.txt")

# Subset columns with normalized count information
DEG_down <- DEG_down[,20:41]

# Make columnes nice for plotting
colnames(DEG_down) <- gsub("\\..*","",colnames(DEG_down))
data <- as.data.frame(colMeans(DEG_down))
colnames(data) <- c("Deg")

# Count H3K27ac in MYC non regulated MYC enhancers
## annotatePeaks.pl MED1_MYC_Enh_No.bed hg38 -d 76F2V.bam.TD A493.bam.TD AU565.bam.TD BT549.bam.TD H128.bam.TD H2171.bam.TD HCC1937.bam.TD HCC1954.bam.TD HCT116.bam.TD HeLa.bam.TD MB231.bam.TD MB361.bam.TD MB436.bam.TD MB468.bam.TD MCF7.bam.TD MM1S.bam.TD SKBR3.bam.TD T47D.bam.TD U2OS.bam.TD U87.bam.TD UACC812.bam.TD ZR751.bam.TD -size 1000 -cpu 24

# Load into R
DEG_no <- read.delim("Cancer_MED1_MYC_Enh_no_count.txt")

# Subset columns with normalized count information
DEG_no <- DEG_no[,20:41]

# Make columnes nice for plotting
colnames(DEG_no) <- gsub("\\..*","",colnames(DEG_no))
data$No <- (colMeans(DEG_no))

# Make function to scale data and add 0 as a starting/ending point
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


# Do the extact same for MYC activated and non regulated promoter regions
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


################# Figure 3e ################# 
# Footprinting analysis was performed accordingly to https://github.com/loosolab/TOBIAS
# Nat Comm. 2020 Bentsen et al.,

# Load in footprint result for ATAC-seq +/- KJ-PYR in BT549 perofmred in MYC regions
MYC_Enh_foot <- read.delim("bindetect_results.txt")

# Add annotation
MYC_Enh_foot$diffexpressed <- "NO"
MYC_Enh_foot$diffexpressed[MYC_Enh_foot$DMSO_KJ_pvalue < 1e-50  & MYC_Enh_foot$DMSO_KJ_change < -0.1] <- "UP"
MYC_Enh_foot$diffexpressed[MYC_Enh_foot$DMSO_KJ_pvalue < 1e-50 & MYC_Enh_foot$DMSO_KJ_change > 0.1] <- "DOWN"

# Make list of factors to highlight in the plot
highlight<-(c("STAT3","FOXC1")) 

# Make the volcano plot
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



################# Figure 3f ################# 
protName<-"MED1" #Set genesymbol for bait protein
protID<-"Q15648" #Set uniport id for bait protein
a<-read.csv("QEHF3_094_RIME19_20210602_TR_peptides_new.csv", sep =";") #Column names for each sample should match the metadata

meta<-read.delim("metafil_new.txt",h=T, stringsAsFactors = F) #Import metadata

n<-11 #Number of samples in TMT experiment
fasta<-"MED1.fasta.txt" #Path to fasta file for bait protein
scaleFunc<-median #Median scaling is not too hard on the data. This is done within groups.
con<-c(MycI_vs_Crtl = "MycI - Crtl") #Describe contrast
padj<-0.2
lfc<-0
highlight<-(c("SOX13","YAP1")) #Proteins to highlight in plots


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
pcaPlot(MSnset_data, labelColumn = "SampleName", pointsize = 2,labelsize = 3) + xlim(-600,200) #PCA plot




#Scale peptide intensities according to group median (DEX/Eth vs IgG)
MSnset_norm_gs <- groupScaling(MSnset_data, 
                               groupingColumn = "Grp")

intensityPlot(MSnset_norm_gs, title = "Peptide intensity distribution")
intensityBoxplot(MSnset_norm_gs, title = "Peptide intensity distribution")
corrPlot(MSnset_norm_gs) 
hierarchicalPlot(MSnset_norm_gs) #Tree plot 
pcaPlot(MSnset_norm_gs, labelColumn = "BioRep", pointsize = 5,labelsize = 0, omitIgG = TRUE, transform = FALSE)




#Retrieve annotation for conversion of peptide info to protein info - NOT normalized data
proteins <- unique(fData(MSnset_norm_gs)$Accessions)
columns <- c("ENTRY-NAME", "PROTEIN-NAMES", "GENES")
data(human_anno)

#Summarise peptide intensities for proteins
MSnset_Pnorm <- summarizeIntensities(MSnset_norm_gs, mean, human_anno)


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

highlight<-(c("STAT3","KDM3A")) 
spec2$highlight <- "NO"
spec2$highlight[spec2$log2FC > 0 & spec2$adj.P.Val < 0.05] <- "YES"
spec2$highlight[spec2$log2FC < 0 & spec2$adj.P.Val < 0.05] <- "YES"


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



################# Figure 3g ################# 
# TOBIAS PlotAggregate --TFBS STAT3_dep_MYC_enh.bed --signals Merged_DMSO_corrected.bw Merged_KJ_corrected.bw --output STAT3_dep_MYC.pdf --share_y both --plot_boundaries --signal-on-x --smooth 10 --flank 120

################# Figure 3h ################# 
# Run on solexa
## annotatePeaks STAT3_dep_MYC_enh.bed hg38 -d combined_MED1_DMSO.TD combined_MED1_KJ.TD -size 5000 -hist 50 -cpu 12

# Load in avarage data
agg <- read.delim("STAT3_dep_MYC_enh_avg.txt", header = TRUE)

# Subset columns for plotting, re-name and melt for plotting
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_prom)<-c("dist","MED1_DMSO","MED1_KJ")
plot <-melt(avr2_prom, id="dist")

# Plot
ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=1)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts")


################# Supp 3d ################## 

# Create dataframe for ERCC-spike in percent for 



################# Supp 4a ################## 
# MYC_activated_enh_motif.txt
# Load in avarage data
Activated <- read.delim("MYC_activated_enh_motif.txt", header = TRUE)
Constant <- read.delim("MYC_non_enh_motif.txt", header = TRUE)

# Subset columns for plotting, re-name and melt for plotting
Activated <-Activated[,c(1,grep("c.Myc_1miss.total.sites", colnames(Activated)))]
colnames(Activated)<-c("dist","Activated")

Constant <-Constant[,c(1,grep("c.Myc_1miss.total.sites", colnames(Constant)))]
colnames(Constant)<-c("dist","Constant")

plot <- cbind(Activated, Constant)
plot <- plot[,c("Activated","Constant","dist")]

plot <-melt(plot, id="dist")

# Plot
ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=1)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + coord_cartesian(ylim=c(0.0008,0.0017)) + xlim(-2400,2400)


################# Supp 4b ################## 
Activated <- read.delim("MED1_MYC_Enh_Up.agg.txt", header = TRUE)
Repressed <- read.delim("MED1_MYC_Enh_Down.agg.txt", header = TRUE)
Constant <- read.delim("MED1_MYC_Enh_No.agg.txt", header = TRUE)

# Subset columns for plotting, re-name and melt for plotting
Activated <-Activated[,c(1,grep("Coverage", colnames(Activated)))]
colnames(Activated)<-c("dist","Activated")

Repressed <-Repressed[,c(1,grep("Coverage", colnames(Repressed)))]
colnames(Repressed)<-c("dist","Repressed")

Constant <-Constant[,c(1,grep("Coverage", colnames(Constant)))]
colnames(Constant)<-c("dist","Constant")

plot <- cbind(Activated, Repressed, Constant)
plot <- plot[,c("Activated","Repressed","Constant","dist")]

plot <-melt(plot, id="dist")

# Plot
ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=1)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts")


################# Supp 4c ################## 
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
################# Supp 4d ################# 
## Giggle analysis BT549 
Giggle_down <- read.csv("Giggle_down.csv")

ggplot(Giggle_down, aes(y=GIGGLE_score, x=reorder(Factor,-GIGGLE_score))) + geom_bar(stat="summary") + geom_point(size = 4) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "NA") + xlab("") + ylab("GIGGLE score") + scale_fill_viridis_d(option="E") + ylim(0,2000)


## Giggle analysis MCF7 
Giggle_down <- read.csv("Giggle_Down_MCF.csv")

ggplot(Giggle_down, aes(y=GIGGLE_score, x=reorder(Factor,-GIGGLE_score)) + geom_bar(stat="summary") + geom_point(size = 4) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "NA") + xlab("") + ylab("GIGGLE score") + scale_fill_viridis_d(option="E") + ylim(0,2000)





################# Supp 4e ################# 
# Load in footprinting score from MYC enhancer and promoter
STAT3_en <- read.delim("STAT3_MA0144.2_enh.txt")
STAT3_pr <- read.delim("STAT3_MA0144.2_prom.txt")

# Make boxplot of the average footprintint score change upon KJ-pyr-9 treatment
ggplot() + 
  geom_boxplot(data=STAT3_en, aes(y=-DMSO_KJ_log2fc, x = "STAT3_en"), outlier.alpha = 0) +
  geom_boxplot(data=STAT3_pr, aes(y=-DMSO_KJ_log2fc, x = "STAT3_pr"), outlier.alpha = 0) +  xlab("") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") +
  geom_vline(xintercept = 0, linetype="dashed") + ylab("Log2Foldchange") + geom_hline(yintercept = 0, linetype="dashed")+ coord_cartesian(ylim=c(-1.5,1.5))

# Perform stastiatical test
t.test(STAT3_en$DMSO_KJ_log2fc, STAT3_pr$DMSO_KJ_log2fc)$p.value

################# Supp 4f ################# 
agg <- read.delim("ER_MYC_enh_avg.txt")

# Grep covarage from regions
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]

# Insert correct names
colnames(avr2_prom)<-c("dist","KJ","DMSO")

# Prepare dataframe for plotting
plot <-melt(avr2_prom, id="dist")

# Plot the average plot
ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=1) +
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts")





