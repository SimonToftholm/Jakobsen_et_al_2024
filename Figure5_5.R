##############################################
################## FIGURE 5 ##################
##############################################

################# Figure 5A ################## 
#### Pause release at enhancer ####
diff_intergenic <- read.delim("MYC_regulated_enhancers_BT549.txt")

enhancer <- read.delim("MED1_MYC_Enh_Down_S2p_pol_count.txt")
enhancer_down <- diff_intergenic[diff_intergenic$diffexpressed == "DOWN",]
enhancer <- merge(enhancer, enhancer_down, by.x = "PeakID..cmd.annotatePeaks.pl.test_down.bed.hg38..d.Combined_S2p_DMOS.TD.Combined_S2p_KJ.TD.Combined_S5p_DMSO.TD.Combined_S5p_KJ.TD....BT549_MED1_KJ.tagdig.combined_MED1_DMSO....BT549_MED1_KJ.tagdig.combined_MED1_KJ..size.2000..cpu.10.",
                  by.y="PeakID..cmd.annotatePeaks.pl.MYC_sites.bed.hg38.")                                 

enhancer$DMSO <-enhancer$Combined_S5p_DMSO.TD.Tag.Count.in.2000.bp..90136786.5.Total..normalization.factor...0.11..effective.total...10000000./enhancer$Combined_S2p_DMOS.TD.Tag.Count.in.2000.bp..92787374.5.Total..normalization.factor...0.11..effective.total...10000000.
enhancer$KJ <- enhancer$Combined_S5p_KJ.TD.Tag.Count.in.2000.bp..66462881.5.Total..normalization.factor...0.15..effective.total...10000000./enhancer$Combined_S2p_KJ.TD.Tag.Count.in.2000.bp..78571007.5.Total..normalization.factor...0.13..effective.total...10000000.

ggplot() + geom_point(data=enhancer, aes(log2(DMSO), y=rank(log2(DMSO))), col = "red") +
          geom_point(data=enhancer, aes(log2(KJ), y=rank(log2(KJ))), col="blue") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")

#### Pause release of MYC promoter regulated genes ####
gene_prom <- read.delim("tss_1kb_hg38_S5P.txt")
gene_body <- read.delim("tss_genebody_hg38_S2P.txt")

gene <- merge(gene_prom, gene_body, by.x = "PeakID..cmd.annotatePeaks.pl.tss.hg38..d.Combined_S2p_DMOS.TD.Combined_S2p_KJ.TD.Combined_S5p_DMSO.TD.Combined_S5p_KJ.TD.New_Combined_S2P_DMSO.TD.New_Combined_S2P_KJ.TD.New_Combined_S5P_DMSO.TD.New_Combined_S5P_KJ.TD....MYC_MED1_ChIP.tagdir.Combined_MYC_BT549...size.1000.", by.y="PeakID..cmd.annotatePeaks.pl.tss.hg38..d.Combined_S2p_DMOS.TD.Combined_S2p_KJ.TD..size.0.10000.")

KJ_Gene <- read.table("DE_analysis_KJ_BT549_DEX_new.txt")
KJ_Gene_MYC_down <- KJ_Gene[KJ_Gene$Log2FC_MYC_vs_Crtl < 0 & KJ_Gene$padj_MYC_vs_Crtl > 0.05,]
gene_KJ_Gene_MYC_down <- merge(KJ_Gene_MYC_down, gene, by.x = "symbol",by.y="Gene.Name.x")
gene_KJ_Gene_MYC_down <- gene_KJ_Gene_MYC_down[!duplicated(gene_KJ_Gene_MYC_down$symbol),]

gene_KJ_Gene_MYC_down <- gene_KJ_Gene_MYC_down[gene_KJ_Gene_MYC_down$...MYC_MED1_ChIP.tagdir.Combined_MYC_BT549..Tag.Count.in.1000.bp..40939315.0.Total..normalization.factor...0.24..effective.total...10000000. > 30,]


gene_KJ_Gene_MYC_down$DMSO <-gene_KJ_Gene_MYC_down$Combined_S5p_DMSO.TD.Tag.Count.in.1000.bp..90136786.5.Total..normalization.factor...0.11..effective.total...10000000./gene_KJ_Gene_MYC_down$Combined_S2p_DMOS.TD.Tag.Count.in.10000.bp..92787374.5.Total..normalization.factor...0.11..effective.total...10000000.
gene_KJ_Gene_MYC_down$KJ <- gene_KJ_Gene_MYC_down$Combined_S5p_KJ.TD.Tag.Count.in.1000.bp..66462881.5.Total..normalization.factor...0.15..effective.total...10000000./gene_KJ_Gene_MYC_down$Combined_S2p_KJ.TD.Tag.Count.in.10000.bp..78571007.5.Total..normalization.factor...0.13..effective.total...10000000.


ggplot() + geom_point(data=gene_KJ_Gene_MYC_down, aes(log2(DMSO), y=rank(DMSO)), col = "red") +
  geom_point(data=gene_KJ_Gene_MYC_down, aes(log2(KJ), y=rank(KJ)), col="blue")  +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")


#### Pause release of MYC regulated genes without promoter binding ####
gene_prom <- read.delim("tss_1kb_hg38_S5P.txt")
gene_body <- read.delim("tss_genebody_hg38_S2P.txt")

gene <- merge(gene_prom, gene_body, by.x = "PeakID..cmd.annotatePeaks.pl.tss.hg38..d.Combined_S2p_DMOS.TD.Combined_S2p_KJ.TD.Combined_S5p_DMSO.TD.Combined_S5p_KJ.TD.New_Combined_S2P_DMSO.TD.New_Combined_S2P_KJ.TD.New_Combined_S5P_DMSO.TD.New_Combined_S5P_KJ.TD....MYC_MED1_ChIP.tagdir.Combined_MYC_BT549...size.1000.", by.y="PeakID..cmd.annotatePeaks.pl.tss.hg38..d.Combined_S2p_DMOS.TD.Combined_S2p_KJ.TD..size.0.10000.")

KJ_Gene <- read.table("DE_analysis_KJ_BT549_DEX_new.txt")
KJ_Gene_MYC_down <- KJ_Gene[KJ_Gene$padj_MYC_vs_Crtl > 0.05,]
gene_KJ_Gene_MYC_down <- merge(KJ_Gene_MYC_down, gene, by.x = "symbol",by.y="Gene.Name.x")
gene_KJ_Gene_MYC_down <- gene_KJ_Gene_MYC_down[!duplicated(gene_KJ_Gene_MYC_down$symbol),]

gene_KJ_Gene_MYC_down <- gene_KJ_Gene_MYC_down[gene_KJ_Gene_MYC_down$...MYC_MED1_ChIP.tagdir.Combined_MYC_BT549..Tag.Count.in.1000.bp..40939315.0.Total..normalization.factor...0.24..effective.total...10000000. < 20,]


gene_KJ_Gene_MYC_down$DMSO <-gene_KJ_Gene_MYC_down$Combined_S5p_DMSO.TD.Tag.Count.in.1000.bp..90136786.5.Total..normalization.factor...0.11..effective.total...10000000./gene_KJ_Gene_MYC_down$Combined_S2p_DMOS.TD.Tag.Count.in.10000.bp..92787374.5.Total..normalization.factor...0.11..effective.total...10000000.
gene_KJ_Gene_MYC_down$KJ <- gene_KJ_Gene_MYC_down$Combined_S5p_KJ.TD.Tag.Count.in.1000.bp..66462881.5.Total..normalization.factor...0.15..effective.total...10000000./gene_KJ_Gene_MYC_down$Combined_S2p_KJ.TD.Tag.Count.in.10000.bp..78571007.5.Total..normalization.factor...0.13..effective.total...10000000.


ggplot() + geom_point(data=gene_KJ_Gene_MYC_down, aes(log2(DMSO), y=rank(DMSO)), col = "red") +
  geom_point(data=gene_KJ_Gene_MYC_down, aes(log2(KJ), y=rank(KJ)), col="blue")


write.table(gene_KJ_Gene_MYC_down[28:30], "MYC_regulated_gene_not_bound.bed", row.names = F, col.names = F, quote = F, sep = "\t")

################# Figure 5B ################## 
#### Plotting average plot over gene sigificant down-regulated with MYC binding ####
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

#### Plotting average plot over enhancer sigificant down-regulated #####
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

################# Figure 5C ################# 
# CAGE-seqquantile normalized average around lost enhancer sites
agg <- read.delim("CAGE_down_count.txt", header = TRUE)
avr2_prom <-agg[,c(1,grep("Tags", colnames(agg)))]
colnames(avr2_prom)<-c("dist","DMSO+","DMSO-","KJ+","KJ-")

library(reshape)
plot <- avr2_prom[,c("dist","DMSO+","DMSO-","KJ+","KJ-")]
plot <-melt(plot, id="dist")

ggplot(plot, aes(x=dist,y=value, col=variable))+
  geom_line(size=1)+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  geom_vline(xintercept = 0, linetype="dashed") + ylab("TagCounts") 




#Rename data packages
library(pheatmap)
library(ggseqlogo)
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
source("CAGEfightR_extensions/enhancers.R")

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









################# Figure 5D################# 
