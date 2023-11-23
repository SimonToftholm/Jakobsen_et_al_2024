##############################################
################## FIGURE 5 ##################
##############################################

# Pre-proccessed data can be downloaded from https://zenodo.org/deposit/8323614

################# Packages #################
library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
################# Figure 5a ################## 
#### Pause release at enhancer
## Run the following on cluster
# annotatePeaks.pl test_down.bed hg38 -d Combined_S2p_DMSO.TD Combined_S2p_KJ.TD Combined_S5p_DMSO.TD Combined_S5p_KJ.TD -size 2000 -cpu 10 > MED1_MYC_Enh_Down_S2p_pol_count.txt

# Load in file of count data
enhancer <- read.delim("MED1_MYC_Enh_Down_S2p_pol_count.txt")
                              
# Calculate ratio between S5p (paused) and S2p (elongating) RNA Pol II
enhancer$DMSO <-enhancer$Combined_S5p_DMSO.TD.Tag.Count.in.2000.bp..90136786.5.Total..normalization.factor...0.11..effective.total...10000000./enhancer$Combined_S2p_DMOS.TD.Tag.Count.in.2000.bp..92787374.5.Total..normalization.factor...0.11..effective.total...10000000.
enhancer$KJ <- enhancer$Combined_S5p_KJ.TD.Tag.Count.in.2000.bp..66462881.5.Total..normalization.factor...0.15..effective.total...10000000./enhancer$Combined_S2p_KJ.TD.Tag.Count.in.2000.bp..78571007.5.Total..normalization.factor...0.13..effective.total...10000000.

# Rank based on ratio and plot
ggplot() + geom_point(data=enhancer, aes(log2(DMSO), y=rank(log2(DMSO))), col = "red") +
          geom_point(data=enhancer, aes(log2(KJ), y=rank(log2(KJ))), col="blue") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")

# Calculate p-value
t.test(enhancer$DMSO,enhancer$KJ)$p.value

#### Pause release of MYC promoter regulated genes
# annotatePeaks.pl tss hg38 -d Combined_S2p_DMOS.TD Combined_S2p_KJ.TD Combined_S5p_DMSO.TD Combined_S5p_KJ.TD New_Combined_S2P_DMSO.TD New_Combined_S2P_KJ.TD New_Combined_S5P_DMSO.TD New_Combined_S5P_KJ.TD ../MYC_MED1_ChIP/tagdir/Combined_MYC_BT549/ -size 1000 > tss_1kb_hg38_S5P.txt
gene_prom <- read.delim("tss_1kb_hg38_S5P.txt")

# annotatePeaks.pl tss hg38 -d Combined_S2p_DMOS.TD Combined_S2p_KJ.TD -size 0,10000 > tss_genebody_hg38_S2P.txt
gene_body <- read.delim("tss_genebody_hg38_S2P.txt")

# Merge the two objects
gene <- merge(gene_prom, gene_body, by.x = "PeakID..cmd.annotatePeaks.pl.tss.hg38..d.Combined_S2p_DMOS.TD.Combined_S2p_KJ.TD.Combined_S5p_DMSO.TD.Combined_S5p_KJ.TD.New_Combined_S2P_DMSO.TD.New_Combined_S2P_KJ.TD.New_Combined_S5P_DMSO.TD.New_Combined_S5P_KJ.TD....MYC_MED1_ChIP.tagdir.Combined_MYC_BT549...size.1000.", by.y="PeakID..cmd.annotatePeaks.pl.tss.hg38..d.Combined_S2p_DMOS.TD.Combined_S2p_KJ.TD..size.0.10000.")

# Add RNA-seq information to get only MYC activated genes
KJ_Gene <- read.table("DE_analysis_KJ_BT549_DEX_new.txt")
KJ_Gene_MYC_down <- KJ_Gene[KJ_Gene$Log2FC_MYC_vs_Crtl < 0 & KJ_Gene$padj_MYC_vs_Crtl < 0.05,]
gene_KJ_Gene_MYC_down <- merge(KJ_Gene_MYC_down, gene, by.x = "symbol",by.y="Gene.Name.x")
gene_KJ_Gene_MYC_down <- gene_KJ_Gene_MYC_down[!duplicated(gene_KJ_Gene_MYC_down$symbol),]

# Take MYC activated genes with high MYC binding at the promoter
gene_KJ_Gene_MYC_down <- gene_KJ_Gene_MYC_down[gene_KJ_Gene_MYC_down$...MYC_MED1_ChIP.tagdir.Combined_MYC_BT549..Tag.Count.in.1000.bp..40939315.0.Total..normalization.factor...0.24..effective.total...10000000. > 30,]

# Calculate ratio of RNA Pol II subtypes
gene_KJ_Gene_MYC_down$DMSO <-gene_KJ_Gene_MYC_down$Combined_S5p_DMSO.TD.Tag.Count.in.1000.bp..90136786.5.Total..normalization.factor...0.11..effective.total...10000000./gene_KJ_Gene_MYC_down$Combined_S2p_DMOS.TD.Tag.Count.in.10000.bp..92787374.5.Total..normalization.factor...0.11..effective.total...10000000.
gene_KJ_Gene_MYC_down$KJ <- gene_KJ_Gene_MYC_down$Combined_S5p_KJ.TD.Tag.Count.in.1000.bp..66462881.5.Total..normalization.factor...0.15..effective.total...10000000./gene_KJ_Gene_MYC_down$Combined_S2p_KJ.TD.Tag.Count.in.10000.bp..78571007.5.Total..normalization.factor...0.13..effective.total...10000000.

# Ranked based on ratio and plot
ggplot() + geom_point(data=gene_KJ_Gene_MYC_down, aes(log2(DMSO), y=rank(DMSO)), col = "red") +
  geom_point(data=gene_KJ_Gene_MYC_down, aes(log2(KJ), y=rank(KJ)), col="green")  +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")

# Calculate p-value
t.test(gene_KJ_Gene_MYC_down$DMSO, gene_KJ_Gene_MYC_down$KJ)$p.value



################# Figure 5b ################## 
#### Plotting average plot over gene sigificant down-regulated with MYC binding

# annotatePeaks.pl tss hg38 -list ../BT549_KJ_POL/MYC_Gene_down.txt -d Combined_S2p_DMSO.test Combined_S2p_KJ.test Combined_S5p_DMSO.test Combined_S5p_KJ.test -size 5000 -hist 50 > test_pol_gene.txt
Enhancer <- read.delim("test_pol_gene.txt")
avr2_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr2_enh)<-c("dist","S2p_D","S2p_K","S5p_D","S5p_K")
avr2_enh <- avr2_enh[,c("dist","S5p_D","S5p_K")]
plot <-melt(avr2_enh, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") 

# annotatePeaks.pl rna hg38 -d Combined_S2p_DMOS.TD Combined_S2p_KJ.TD Combined_S5p_DMSO.TD Combined_S5p_KJ.TD -hist 100 -cpu 10 -list MYC_Gene_down.txt > POL_gene_body_agg.txt
Enhancer <- read.delim("POL_gene_body_agg.txt")
avr2_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr2_enh)<-c("dist","S2p_D","S2p_K","S5p_D","S5p_K")
avr2_enh <- avr2_enh[,c("dist","S2p_D","S2p_K")]
plot <-melt(avr2_enh, id="dist") 

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + ylim(0,6)

#### Plotting average plot over enhancer sigificant down-regulated
# annotatePeaks.pl ../BT549_MED1_KJ/MED1_MYC_Enh_Down.bed hg38 -d Combined_S2p_DMSO.test Combined_S2p_KJ.test Combined_S5p_DMSO.test Combined_S5p_KJ.test -size 5000 -hist 50 > test_pol.txt
Enhancer <- read.delim("test_pol.txt")
avr2_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr2_enh) <- c("dist","S2p_D","S2p_K","S5p_D","S5p_K")
avr2_enh <- avr2_enh[,c("dist","S2p_D","S2p_K")]
plot <-melt(avr2_enh, id="dist") 

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + ylim(0.5,2.5)

# annotatePeaks.pl ../BT549_MED1_KJ/MED1_MYC_Enh_Down.bed hg38 -d Combined_S2p_DMSO.test Combined_S2p_KJ.test Combined_S5p_DMSO.test Combined_S5p_KJ.test -size 5000 -hist 50 > test_pol.txt
Enhancer <- read.delim("test_pol.txt")
avr2_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr2_enh) <- c("dist","S2p_D","S2p_K","S5p_D","S5p_K")
avr2_enh <- avr2_enh[,c("dist","S5p_D","S5p_K")]
plot <-melt(avr2_enh, id="dist") 

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2)+
  scale_fill_viridis_d(option="E") + scale_color_viridis_d(option="E") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + ylim(0.5,2.5)

################# Figure 5c ################# 
# Load diffbind analyse from RNAPII +/- the three inhibitors

#Read in the samplesheet (Run on server)
#samples<-read.csv("SampleSheetDiffBind_POL_NewInh.csv", sep = ";")
#samples <- samples[1:11,1:10]
#bt<-dba(sampleSheet = samples)

#bt.cons<- dba.peakset(bt, consensus=DBA_CONDITION, minOverlap=0.99)
#bt.cons <- dba(bt.cons, mask=bt.cons$masks$Consensus, minOverlap=1)
# dba.plotMA(bt.cons.count, contrast = 1, method=DBA_DESEQ2, bUsePval = 0.05)

#bt.cons.count<-dba.count(bt, peaks = consensus_peaks, bParallel = TRUE, filter = 0, summits = FALSE )

# Create datafrom containing the information from the DESeq2 analysis
#bt.cons.count <- dba.contrast(bt.cons.count, categories=DBA_CONDITION, minMembers = 2)
#bt.cons.count <- dba.analyze(bt.cons.count, method=DBA_DESEQ2, bBlacklist=F, bGreylist=F)

load("RNAPII_MYC_regions.image.Rdata")

KJ_report<-dba.report(bt.cons.count, th=1, DataType = DBA_DATA_FRAME, method = DBA_DESEQ2, contrast = 10)
KJ_rpkm<-bt.cons.count$peaks[[1]] #Get additional count info for each peak
KJ_rpkm$id<-rownames(KJ_rpkm)
KJ_report$id<-rownames(KJ_report)
KJ_data<-merge(KJ_report, KJ_rpkm, by="id")
colnames(KJ_data)[2:ncol(KJ_data)]<-paste("KJ", colnames(KJ_data)[2:ncol(KJ_data)],sep = "_") #Add KJ to all column names

MYCMI6_report<-dba.report(bt.cons.count, th=1, DataType = DBA_DATA_FRAME, method = DBA_DESEQ2, contrast = 1)
MYCMI6_rpkm<-bt.cons.count$peaks[[1]] #Get additional count info for each peak
MYCMI6_rpkm$id<-rownames(MYCMI6_rpkm)
MYCMI6_report$id<-rownames(MYCMI6_report)
MYCMI6_data<-merge(MYCMI6_report, MYCMI6_rpkm, by="id")
colnames(MYCMI6_data)[2:ncol(MYCMI6_data)]<-paste("MYCMI6", colnames(MYCMI6_data)[2:ncol(MYCMI6_data)],sep = "_") #Add MYCMI6 to all column names

MYCi975_report<-dba.report(bt.cons.count, th=1, DataType = DBA_DATA_FRAME, method = DBA_DESEQ2, contrast = 2)
MYCi975_rpkm<-bt.cons.count$peaks[[1]] #Get additional count info for each peak
MYCi975_rpkm$id<-rownames(MYCi975_rpkm)
MYCi975_report$id<-rownames(MYCi975_report)
MYCi975_data<-merge(MYCi975_report, MYCi975_rpkm, by="id")
colnames(MYCi975_data)[2:ncol(MYCi975_data)]<-paste("MYCi975", colnames(MYCi975_data)[2:ncol(MYCi975_data)],sep = "_") #Add MYCi975 to all column names

# Load diffbind analyse from MED1 +/- KJ
load("MED1_KJ_MYC_regions.image.Rdata")

# Set contrast and run DEseq2 analysis
bt.cons.count <- dba.contrast(bt.cons.count, categories=DBA_CONDITION, minMembers = 2)
bt.cons.count <- dba.analyze(bt.cons.count, method=DBA_DESEQ2, bBlacklist=F, bGreylist=F)

# Create datafrom containing the information from the DESeq2 analysis
MYC_report<-dba.report(bt.cons.count, th=1, DataType = DBA_DATA_FRAME, method = DBA_DESEQ2, contrast = 1, bAll = TRUE)
MYC_rpkm<-bt.cons.count$peaks[[1]] #Get additional count info for each peak
MYC_rpkm$id<-rownames(MYC_rpkm)
MYC_report$id<-rownames(MYC_report)
MYC_data<-merge(MYC_report, MYC_rpkm, by="id")
colnames(MYC_data)[2:ncol(MYC_data)]<-paste("MYC", colnames(MYC_data)[2:ncol(MYC_data)],sep = "_") #Add MYC to all column names


merged <- merge(KJ_data, MYC_data, by = "id")
merged <- merge(MYCMI6_data, merged, by = "id")
merged <- merge(MYCi975_data, merged, by = "id")

ann <- read.delim("MYC_sites.ann.txt")
MYC_data_ann <- merge(ann, merged, by.x ="PeakID..cmd.annotatePeaks.pl.MYC_sites.bed.hg38.", by.y="id")

### Volcano plot for enhancer and promoter regions BT549
MYC_data_ann$Ann <- sub(" .*","",MYC_data_ann$Annotation)
diff_intergenic <- MYC_data_ann[MYC_data_ann$Ann == "Intergenic" | MYC_data_ann$Ann == "intron",]
diff_promoter <- MYC_data_ann[MYC_data_ann$Ann == "promoter-TSS",]

# Sig enhancer 
diff_intergenic$diffexpressed <- "NO"
diff_intergenic$diffexpressed[diff_intergenic$KJ_Fold < 0 & diff_intergenic$`KJ_p-value` < 0.05] <- "KJ_lost"
diff_intergenic$diffexpressed[diff_intergenic$KJ_Fold > 0 & diff_intergenic$`KJ_p-value` < 0.05] <- "KJ_gained"
diff_intergenic_sig <- diff_intergenic[!diff_intergenic$diffexpressed == "NO",]

ggplot() + 
  geom_point(data=diff_intergenic, aes(x=MYC_Fold, y=-log(MYC_FDR)), color = "grey", alpha = 1) +
  geom_point(data=diff_intergenic_sig, aes(x=MYC_Fold, y=-log(MYC_FDR), col = diffexpressed), alpha = 1) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("-log10(padj)") + xlab("log2(KJ/DMSO)") + ylim(0,50) + xlim(-10,10) + scale_color_manual(values = c("#B11518","#253884"))

table(diff_intergenic$diffexpressed)


# Sig enhancer 
diff_intergenic$diffexpressed <- "NO"
diff_intergenic$diffexpressed[diff_intergenic$MYCMI6_Fold > 0 & diff_intergenic$`MYCMI6_p-value` < 0.05] <- "MYCMI6_lost"
diff_intergenic$diffexpressed[diff_intergenic$MYCMI6_Fold < 0 & diff_intergenic$`MYCMI6_p-value` < 0.05] <- "MYCMI6_gained"
diff_intergenic_sig <- diff_intergenic[!diff_intergenic$diffexpressed == "NO",]

ggplot() + 
  geom_point(data=diff_intergenic, aes(x=MYC_Fold, y=-log(MYC_FDR)), color = "grey", alpha = 1) +
  geom_point(data=diff_intergenic_sig, aes(x=MYC_Fold, y=-log(MYC_FDR), col = diffexpressed), alpha = 1) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("-log10(padj)") + xlab("log2(KJ/DMSO)") + ylim(0,50) + xlim(-10,10) + scale_color_manual(values = c("#B11518","#253884"))

table(diff_intergenic$diffexpressed)

# Sig enhancer 
diff_intergenic$diffexpressed <- "NO"
diff_intergenic$diffexpressed[diff_intergenic$MYCi975_Fold > 0 & diff_intergenic$`MYCi975_p-value` < 0.05] <- "MYCi975_lost"
diff_intergenic$diffexpressed[diff_intergenic$MYCi975_Fold < 0 & diff_intergenic$`MYCi975_p-value` < 0.05] <- "MYCi975_gained"
diff_intergenic_sig <- diff_intergenic[!diff_intergenic$diffexpressed == "NO",]

ggplot() + 
  geom_point(data=diff_intergenic, aes(x=MYC_Fold, y=-log(MYC_FDR)), color = "grey", alpha = 1) +
  geom_point(data=diff_intergenic_sig, aes(x=MYC_Fold, y=-log(MYC_FDR), col = diffexpressed), alpha = 1) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + ylab("-log10(padj)") + xlab("log2(KJ/DMSO)") + ylim(0,50) + xlim(-10,10) + scale_color_manual(values = c("#B11518","#253884"))

table(diff_intergenic$diffexpressed)







################# Figure 5d ################# 
agg <- read.delim("U2OS_MYC_enh_hist.txt", header = TRUE)
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_prom)<-c("dist","Dox_2","No2","Dox_5","No_5")

library(reshape)
plot <-melt(avr2_prom, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=1)+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") 


# Total POL II data??

################# Figure 5e ################# 
# CAGE-seqquantile normalized average around lost enhancer sites
agg <- read.delim("CAGE_down_count.txt", header = TRUE)
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_prom)<-c("dist","DMSO","KJ")

library(reshape)
plot <-melt(avr2_prom, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=1)+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlim(-400,400)







################# Supp 6a ################# 
gene_prom <- read.delim("tss_1kb_hg38_S5P.txt")
gene_body <- read.delim("tss_genebody_hg38_S2P.txt")

gene <- merge(gene_prom, gene_body, by.x = "PeakID..cmd.annotatePeaks.pl.tss.hg38..d.Combined_S2p_DMOS.TD.Combined_S2p_KJ.TD.Combined_S5p_DMSO.TD.Combined_S5p_KJ.TD.New_Combined_S2P_DMSO.TD.New_Combined_S2P_KJ.TD.New_Combined_S5P_DMSO.TD.New_Combined_S5P_KJ.TD....MYC_MED1_ChIP.tagdir.Combined_MYC_BT549...size.1000.", by.y="PeakID..cmd.annotatePeaks.pl.tss.hg38..d.Combined_S2p_DMOS.TD.Combined_S2p_KJ.TD..size.0.10000.")

KJ_Gene <- read.table("DE_analysis_KJ_BT549_DEX_new.txt")
KJ_Gene_MYC_down <- KJ_Gene[KJ_Gene$padj_MYC_vs_Crtl > 0.2,]
gene_KJ_Gene_MYC_down <- merge(KJ_Gene_MYC_down, gene, by.x = "symbol",by.y="Gene.Name.x")
gene_KJ_Gene_MYC_down <- gene_KJ_Gene_MYC_down[!duplicated(gene_KJ_Gene_MYC_down$symbol),]

gene_KJ_Gene_MYC_down <- gene_KJ_Gene_MYC_down[gene_KJ_Gene_MYC_down$...MYC_MED1_ChIP.tagdir.Combined_MYC_BT549..Tag.Count.in.1000.bp..40939315.0.Total..normalization.factor...0.24..effective.total...10000000. < 30,]


gene_KJ_Gene_MYC_down$DMSO <-gene_KJ_Gene_MYC_down$Combined_S5p_DMSO.TD.Tag.Count.in.1000.bp..90136786.5.Total..normalization.factor...0.11..effective.total...10000000./gene_KJ_Gene_MYC_down$Combined_S2p_DMOS.TD.Tag.Count.in.10000.bp..92787374.5.Total..normalization.factor...0.11..effective.total...10000000.
gene_KJ_Gene_MYC_down$KJ <- gene_KJ_Gene_MYC_down$Combined_S5p_KJ.TD.Tag.Count.in.1000.bp..66462881.5.Total..normalization.factor...0.15..effective.total...10000000./gene_KJ_Gene_MYC_down$Combined_S2p_KJ.TD.Tag.Count.in.10000.bp..78571007.5.Total..normalization.factor...0.13..effective.total...10000000.


ggplot() + geom_point(data=gene_KJ_Gene_MYC_down, aes(log2(DMSO), y=rank(DMSO)), col = "#253884") +
  geom_point(data=gene_KJ_Gene_MYC_down, aes(log2(KJ), y=rank(KJ)), col="#929BC1")  +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")

# Calculate p-value
t.test(gene_KJ_Gene_MYC_down$DMSO, gene_KJ_Gene_MYC_down$KJ)


################# Supp 6b ################# 
Enhancer <- read.delim("New_inh_MYC_enh_agg.txt")

avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","VEH","MYCMI6","MYCi975")

plot <-melt(avr1_enh, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp") 

Enhancer <- read.delim("New_inh_MYC_prom_agg.txt")

avr1_enh <-Enhancer[,c(1,grep("Coverage", colnames(Enhancer)))]
colnames(avr1_enh)<-c("dist","VEH","MYCMI6","MYCi975")

plot <-melt(avr1_enh, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=2) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") + xlab("bp") 

################# Supp 6c ################# 
agg <- read.delim("MYC_down_Intergenic_ann.txt", header = TRUE)
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


agg <- read.delim("MYC_down_Intron_ann.txt", header = TRUE)
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


################# Supp 6d ################# 


#Rename data packages
library(pheatmap)
library(viridis)
library(magrittr)
library(ggforce)
library(ggthemes)
library(tidyverse)

library(CAGEfightR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(org.Hs.eg.db)
library(GenomicRanges)
library(SummarizedExperiment)
library(GenomicFeatures)
library(BiocParallel)
library(InteractionSet)
library(Gviz)
library(ensembldb)

library(DESeq2)
library(limma)
library(edgeR)
library(statmod)
library(BiasedUrn)
library(sva)

library(TFBSTools)
library(motifmatchr)
library(pathview)
library(EnsDb.Hsapiens.v86)

bsg <-TxDb.Hsapiens.UCSC.hg38.knownGene
TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
odb <- org.Hs.eg.db
genomeInfo <- SeqinfoForUCSCGenome("hg38")
source("enhancers.R")

genomeInfo <- keepStandardChromosomes(SeqinfoForUCSCGenome("hg38"),species="Homo_sapiens")

seqlevels(TxDb) <- c(seqlevels(genomeInfo))  #Keep only standard chromosomes


myc_up <- read.table("MED1_MYC_Enh_Up.bed",header=T)
med_down <- read.table("MED1_MYC_Enh_Down.bed",header=T)
med_un <- read.table("MED1_MYC_Enh_No.bed",header=T)
myc_up$patt <- "up"
med_down$patt <- "down"
med_un$patt <- "no"

rang <- GRanges(rbind(myc_up,med_down,med_un),seqinfo = genomeInfo)

# Make metadata
Name <- c("DMSO1","DMSO_2","DMSO_3","KJ_1","KJ_2","KJ_3")
Class <- c("DMSO","DMSO","DMSO","KJ","KJ","KJ")
Batch <- c("A","B","C","A","B","C")
BigWigPlus <- c("BT549_DMSO_Rep1.plus.bw","BT549_DMSO_Rep2.plus.bw","BT549_DMSO_Rep3.plus.bw","BT549_KJ_Rep1.plus.bw","BT549_KJ_Rep2.plus.bw","BT549_KJ_Rep3.plus.bw")
BigWigMinus <- c("BT549_DMSO_Rep1.minus.bw","BT549_DMSO_Rep2.minus.bw","BT549_DMSO_Rep3.minus.bw","BT549_KJ_Rep1.minus.bw","BT549_KJ_Rep2.minus.bw","BT549_KJ_Rep3.minus.bw")


sampleTable <- data.frame(Class ,Name,Batch, BigWigPlus, BigWigMinus, row.names = Name)
sampleTable

# Load in path for bigwig files
bw_plus <- system.file("extdata", sampleTable$BigWigPlus, package = "CAGEfightR")
bw_minus <- system.file("extdata", sampleTable$BigWigMinus,package = "CAGEfightR")


# Create two named BigWigFileList-objects:
bw_plus <- BigWigFileList(bw_plus)
bw_minus <- BigWigFileList(bw_minus)
names(bw_plus) <- sampleTable$Name
names(bw_minus) <- sampleTable$Name


#Quantify CTSSs
register(MulticoreParam(workers=12))


CTSSs <- quantifyCTSSs(plusStrand = bw_plus,
                       minusStrand = bw_minus,
                       genome=genomeInfo,
                       design = sampleTable)


med <- quantifyClusters(CTSSs,rang)


med <- med %>% calcTPM()

assay(med,"q_TPM") <- normalizeQuantiles(assay(med,"TPM"))


KJ <- calcPooled(subset(med, select = colData(med)$Class=="KJ"),outputColumn="KJ_score",inputAssay="q_TPM")
DMSO <- calcPooled(subset(med, select = colData(med)$Class=="DMSO"),outputColumn="DMSO",inputAssay="q_TPM")

rowData(med)$KJ_score <- rowData(KJ)$KJ_score/3
rowData(med)$DMSO_score <- rowData(DMSO)$DMSO/3

rowData(med)$FC <- log2((rowData(med)$KJ_score)/(rowData(med)$DMSO_score))

#Try plotting
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R")

med %>% rowData %>% as.data.frame %>%
  ggplot(aes(patt,FC,fill=patt)) + geom_flat_violin() + xlab("") + geom_boxplot() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") +
  geom_vline(xintercept = 0, linetype="dashed") + ylab("Log2Foldchange") + coord_flip() + geom_hline(yintercept = 0, linetype="dashed") + ylim(-6,6)





################# Supp 6e ################# 
# annotatePeaks.pl eRNA.bed hg38 -d Ribo_DMSO.TD Ribo_KJ.TD ../../BT549_MED1_KJ/tagdig/combined_MED1_DMSO ../../BT549_MED1_KJ/tagdig/combined_MED1_KJ ../../BT549_ATAC_KJ/combined_ATAC_DMSO ../../BT549_ATAC_KJ/combined_ATAC_KJ ../../CAGEseq/combined_CAGE_DMSO.TD ../../CAGEseq/combined_CAGE_KJ.TD -size 1000 -cpu 12 > eRNA_count.txt

# Load in data 
eRNA_count <- read.delim("eRNA_count.txt")

# Make boxplot
boxplot(eRNA_count$Ribo_DMSO.TD.Tag.Count.in.1000.bp..93528544.0.Total..normalization.factor...0.11..effective.total...10000000., eRNA_count$Ribo_KJ.TD.Tag.Count.in.1000.bp..38739617.5.Total..normalization.factor...0.26..effective.total...10000000., outline = F, ylim = c(0,4))

# Calculate p-value
t.test(eRNA_count$Ribo_DMSO.TD.Tag.Count.in.1000.bp..93528544.0.Total..normalization.factor...0.11..effective.total...10000000., eRNA_count$Ribo_KJ.TD.Tag.Count.in.1000.bp..38739617.5.Total..normalization.factor...0.26..effective.total...10000000., paired = T)

################# Supp 6f ################# 
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
















