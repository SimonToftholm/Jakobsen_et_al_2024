##############################################
################## FIGURE 4 ##################
##############################################

################# Figure 4A ################## 

protName<-"MYC" #Set genesymbol for bait protein
protID<-"P01106" #Set uniport id for bait protein
a<-read.delim("QEHF3_086_RIME16_frak1to9_220920_peptides.txt") #Column names for each sample should match the metadata
a<- a[,1:17]
a <- as.numeric(a)
meta<-read.delim("metafil.txt",h=T, stringsAsFactors = F) #Import metadata
n<-11 #Number of samples in TMT experiment
fasta<-"P01106_MYC_FASTA.txt" #Path to fasta file for bait protein
scaleFunc<-median #Median scaling is not too hard on the data. This is done within groups.
con<-c(DEX_vs_Eth = "Eth - IgG") #Describe contrast
padj<-0.05
lfc<-0.3
highlight<-(c("KAT2A")) #Proteins to highlight in plots


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
library(UniProt.ws)
library(dplyr)
proteins <- unique(fData(MSnset_norm_gs)$Accessions)
columns <- c("ENTRY-NAME", "PROTEIN-NAMES", "GENES")
keyt <- UniProt.ws::UniProt.ws(taxId = 9606)
hs_anno <- UniProt.ws::select(hs, proteins, columns, "UniProtKB") %>%
  as_tibble() %>%
  mutate(GeneSymbol = gsub(" .*", "", GENES)) %>%
  dplyr::select(
    Accessions = "UNIPROTKB", Gene = "ENTRY-NAME",
    Description = "PROTEIN-NAMES", GeneSymbol
  )

#Summarise peptide intensities for proteins
MSnset_Pnorm <- summarizeIntensities(MSnset_norm_gs, sum, hs_anno)

#Get diff results
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




highlight<-(c("KAT2A","STAT3","EP300","MED1","MAX","EHMT2")) 
spec2$highlight <- "NO"



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







#Finding transcriptional regulators in the MYC interactome +/- DEX (BT549)
transReg<-read.delim("Transcriptional_regulators_full_list_from_CROPseq_design_script.txt",h=T)
MYC_transreg<-spec[spec$GeneSymbol %in% transReg$symbol,] # = 591 proteins out of 2689 are transcriptional regulators 
MYC_Dex_transreg_up_sig<-MYC_transreg[MYC_transreg$adj.P.Val<0.05 & MYC_transreg$log2FC>0.3, ] # = 2 protens: NR3C1 and ETV6
MYC_Dex_transreg_up<-MYC_transreg[MYC_transreg$log2FC>0, ]   # = 415 proteins
MYC_Dex_transreg_down<-MYC_transreg[MYC_transreg$log2FC<0, ] # = 165 proteins 


#MYC interactome +/-DEX 
colnames(MYC_transreg)
head(MYC_transreg)
RM_38<-MYC_transreg[,c(4,17,23)]
head(RM_38)



###sammenligning af MYC interactome RM-38 og GR interactome RM-30


# Make vectors of RM_38 and RM_30 transcriptional regulators-lists
v_RM_38_transreg <- as.vector(RM_38$GeneSymbol)
v_RM_30_dex_transreg <- as.vector(GR_dex_transreg$GeneSymbol)
v_RM_29_dex_transreg <- as.vector(RM_29_dex_transreg$GeneSymbol)
v_RM_19_FOXK1_transreg <- as.vector(RM_19_FOXK1_transreg$GeneSymbol)

#install.packages('eulerr') #only do this once 
library(eulerr)

#Show dublicates in DEP vector, should be =0, If not go back and do AAA 
v_RM_38_transreg[duplicated(v_RM_38_transreg)]
v_RM_30_dex_transreg[duplicated(v_RM_30_dex_transreg)]
v_RM_29_dex_transreg[duplicated(v_RM_29_dex_transreg)]
v_RM_19_FOXK1_transreg[duplicated(v_RM_19_FOXK1_transreg)]

# AAA # Remove dublicated lines, problem is, that the analysis removes both 
GR_dex_transreg <- GR_dex_transreg[!GR_dex_transreg$GeneSymbol == "CDKN2A", ]
RM_19_FOXK1_transreg <- RM_19_FOXK1_transreg[!RM_19_FOXK1_transreg$GeneSymbol == "CDKN2A",]

Euler <- eulerr::euler(list("MYC (+/-DEX)" = v_RM_38_transreg, "GR (GR/IgG)" = v_RM_30_dex_transreg, "GR (+/-DEX)" = v_RM_29_dex_transreg, "FOXK1 (-DEX)" = v_RM_19_FOXK1_transreg))

plot(Euler,
     edges = FALSE,
     legend = FALSE,
     quantities = TRUE,
     strips = NULL,
     main = "Venn Diagram - MYC and GR interactome overlap",
     n = 200L,
     adjust_labels = TRUE,
     fills = list(fill = c("blue","steelblue4", "red", "green"), alpha = 0.5),
     labels = list(col = "white", font = 2))

RM38_RM30_overlap<- RM_38[RM_38$GeneSymbol %in% GR_dex_transreg$GeneSymbol,]
RM38_RM30_RM29_overlap<-RM38_RM30_overlap[RM38_RM30_overlap$GeneSymbol %in% RM_29_dex_transreg$GeneSymbol,]




################# Figure 4B ################## 


load("GR.02022022.image.Rdata")
library(DiffBind)

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

shared<-combi[combi$GR_FDR<0.1 & combi$GR_Fold < 0 & combi$GR_Conc_GR_DEX > 5 & combi$MYC_Conc_MYC_DEX > 5,] # Sites were MYC is already precent but GR is recruited
gained_MYC_GR<-shared[shared$`MYC_p-value`<0.1 & shared$MYC_Fold > 0  & shared$GR_Fold < 0,] # GR gained sites where MYC is also recruited
uncahnged_MYC_GR<-shared[shared$`MYC_p-value`>0.1 & shared$GR_Fold < 0,] # GR gained sites where MYC is also recruited


write.table(gained_MYC_GR[,2:4],"Gained_MYC_Diff_DiffBindsites.bed", quote=F, append=F, row.names = F, col.names = F, sep = "\t")
write.table(uncahnged_MYC_GR[,2:4],"Constand_MYC_Diff_DiffBindsites.bed", quote=F, append=F, row.names = F, col.names = F, sep = "\t")



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






################# Figure 4C ################## 
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


#Average plot
agg <- read.delim("Gained_MYC_Diff_DiffBindsites_Enhancer_KJ_MED1_agg.txt", header = TRUE)
avr2_prom <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_prom)<-c("dist","DEX_KJ","DEX","KJ","VEH")

library(reshape)
plot <- avr2_prom[,c("dist","VEH","DEX_KJ","DEX","KJ")]
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


## Same for promoters
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

a<-read.delim("Gained_MYC_Diff_DiffBindsites_promoter_MED1_KJ_prom.txt",h=T) 
a2<-a[order(a$X0.2, decreasing = T),2:ncol(a)] 
MED1_DEX_KJ <-a2[,1:101]
MED1_DEX <-a2[,102:202]
MYC_KJ <-a2[,203:303]
MYC_VEH <-a2[,304:404]

MED1_DEX_KJ_hist <- Heatmap(as.matrix(MED1_DEX_KJ), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,12,15),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MED1_DEX_KJ", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MED1_DEX_hist <- Heatmap(as.matrix(MED1_DEX), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,12,15),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MED1_DEX", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MYC_KJ_hist <- Heatmap(as.matrix(MYC_KJ), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,12,15),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MYC_KJ", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))
MYC_VEH_hist <- Heatmap(as.matrix(MYC_VEH), cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F, col = colorRamp2(c(0,12,15),c("white","blue1","blue4")), use_raster = T, raster_quality = 10, column_split = c(rep("MYC_VEH", 101)), border = T,heatmap_legend_param = list(direction="horizontal", title_position="topcenter"),height = unit(6,"cm"), width = unit(1.5,"cm"))

draw(MYC_VEH_hist+MYC_KJ_hist+MED1_DEX_hist+MED1_DEX_KJ_hist, heatmap_legend_side="bottom",row_title=paste(nrow(a2),"sites", collapse = " "))



################# Figure 4D ################## 

# Promoter enhancer and gene expression
promoter <- read.delim("tss_hg38_MYC_DEX_count.txt")
enhancer <- read.delim("GR_MYC_sites_DESeq_Intergenic.count.txt")

merged <- merge(enhancer, promoter, by = "Gene.Name")

merged <- merged[!duplicated(merged$Start.x),]
KJ_Gene <- read.table("DE_analysis_KJ_BT549_DEX_new.txt")
merged <- merge(merged, KJ_Gene, by.x = "Gene.Name", by.y = "symbol")


MYC_prom_low <- merged[merged$...tagdig.MYC_DEX_merged..Tag.Count.in.1000.bp..49362524.0.Total..normalization.factor...0.20..effective.total...10000000..y < 15, ]
MYC_prom_high <- merged[merged$...tagdig.MYC_DEX_merged..Tag.Count.in.1000.bp..49362524.0.Total..normalization.factor...0.20..effective.total...10000000..y > 15, ]

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

################# Figure 4E ################## 
# Genes expression
MYC_prom_low$type <- rep("low")
MYC_prom_high$type <- rep("high")

plot <- rbind(MYC_prom_low, MYC_prom_high)
plot <- plot[,c("type","Crtl_DEX1")]
plot <- melt(plot, id = "type")


ggplot(plot, aes(x=type, y=(value), fill = type))  + geom_boxplot(alpha = 1, outlier.alpha = 0) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  xlab(" ") + ylab("Normalized counts") + scale_fill_viridis_d(option="E") + coord_cartesian(ylim=c(0,5000))



# Change +/- KJ with DEX
merged_gene_sig <- merged[merged$padj_Crtl.DEX_vs_Crtl < 1,]
merged_gene_sig <- merged_gene_sig[!duplicated(merged_gene_sig$Gene.Name),]
merged_gene_sig <- merged_gene_sig[complete.cases(merged_gene_sig$Log2FC_Crtl.DEX_vs_Crtl),]
MYC_prom_low_sig <- merged_gene_sig[merged_gene_sig$...tagdig.MYC_DEX_merged..Tag.Count.in.1000.bp..49362524.0.Total..normalization.factor...0.20..effective.total...10000000..y < 15 & merged_gene_sig$Log2FC_Crtl.DEX_vs_Crtl > 0 & abs(merged_gene_sig$Distance.to.TSS.x) < 50000, ]
MYC_prom_high_sig <- merged_gene_sig[merged_gene_sig$...tagdig.MYC_DEX_merged..Tag.Count.in.1000.bp..49362524.0.Total..normalization.factor...0.20..effective.total...10000000..y > 15 & merged_gene_sig$Log2FC_Crtl.DEX_vs_Crtl > 0 & abs(merged_gene_sig$Distance.to.TSS.x) < 50000,]


MYC_prom_low_sig$type <- rep("low")
MYC_prom_high_sig$type <- rep("high")

plot <- rbind(MYC_prom_low_sig, MYC_prom_high_sig)
plot <- plot[,c("type","Log2FC_Crtl.DEX_vs_Crtl","Log2FC_MYC_DEX_vs_MYC")]
plot <- melt(plot, id = "type")

ggplot(plot, aes(x=type, y=value, fill = variable))  + geom_boxplot(alpha = 1, outlier.alpha = 0) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  xlab(" ") + ylab("Mean Expression log Intensity") + scale_fill_viridis_d(option="E") + coord_cartesian(ylim=c(-2,3)) + geom_hline(yintercept=0, linetype="dashed", color = "black")

t.test(MYC_prom_high_sig$Log2FC_Crtl.DEX_vs_Crtl,
       MYC_prom_high_sig$Log2FC_MYC_DEX_vs_MYC, paired = T)








heatmap <- rbind(MYC_prom_low_sig, MYC_prom_high_sig)
MYC_prom_high_sig_heat <- heatmap[,c("Crtl1","Crtl2","Crtl3","MYC1","MYC2","MYC3","Crtl_DEX1","Crtl_DEX2","Crtl_DEX3","MYC_DEX1","MYC_DEX2","MYC_DEX3")]
rownames(MYC_prom_high_sig_heat) <- heatmap$Gene.Name


library("gplots")
heatmap.2(as.matrix(MYC_prom_high_sig_heat), scale = "row", col = bluered(100), trace = "none", density.info = "none", Colv=FALSE, Rowv=TRUE)


################# Figure 4F ################## 
rotName<-"MED1" #Set genesymbol for bait protein
protID<-"Q15648" #Set uniport id for bait protein
a<-read.delim("QEHF3_094_RIME19_20210602_TR_peptides.txt",h=T) #Column names for each sample should match the metadata
a<- a[,1:17]
meta<-read.delim("metafil.txt",h=T, stringsAsFactors = F) #Import metadata
n<-10 #Number of samples in TMT experiment
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
spec<-diffexp[diffexp$controlLogFoldChange>0.5,]
high<-spec[spec$GeneSymbol %in% highlight,c("Description","GeneSymbol")]
colnames(high)[2]<-"text"
spec2<-merge(spec,high, by="Description", all.x = T)
spec2$text[is.na(spec2$text)]<-""


highlight<-(c("NR3C1"))
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

################# Figure 4G
## Compare enhancer-genes with TCGA
library(cBioPortalData)
library(ggplot2)
library(survminer)
library(survival)
brca_tcga<-cBioDataPack("brca_tcga", ask = FALSE) #Import everything for TCGA breast cancer
tcga_exp<-brca_tcga[["RNA_Seq_v2_expression_median"]] #Extract mRNA data from TCGA - Firehose legacy data containing RNAseq from 1100 patients. Info on types of data can be found her https://www.cbioportal.org/faq#rna. The downloaded data here is RSEM format.
tcga_exp2<-assay(tcga_exp) #Get expression matrix
tcga_exp2 <-as.data.frame(tcga_exp2)
tcga_clin<-data.frame(colData(brca_tcga)) #Get clinical information

ER_neg <- tcga_clin[tcga_clin$ER_STATUS_BY_IHC == "Negative",] # Get ER postive and negative
ER_pos <- tcga_clin[tcga_clin$ER_STATUS_BY_IHC == "Positive",]

ER_neg_exp <- as.data.frame(tcga_exp2[,colnames(tcga_exp2) %in% ER_neg$SAMPLE_ID]) # Merge with expression matrix
ER_pos_exp <- as.data.frame(tcga_exp2[,colnames(tcga_exp2) %in% ER_pos$SAMPLE_ID])


ER_neg_exp_DEX <- as.data.frame(tcga_exp2[rownames(tcga_exp2) %in% MYC_prom_high_sig$Gene.Name,])

mean <- (colMeans(ER_neg_exp_DEX))
GR <- (tcga_exp2[rownames(tcga_exp2) %in% c("NR3C1"),])
MYC <- (tcga_exp2[rownames(tcga_exp2) %in% c("MYC"),])

mean <- as.data.frame(t(rbind(mean,MYC,GR)))

ggplot(mean, aes(x=rank(MYC), y=mean)) + geom_smooth() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  xlab("Ranked MYC expression") + ylab("Mean Expression")




agg <- read.delim("Shared_DiffBindsites_agg.txt", header = TRUE)
avr2_yes <-agg[,c(1,grep("Coverage", colnames(agg)))]
colnames(avr2_yes)<-c("dist","GR_DEX","GR","MED1_DEX","MED1","MYC_DEX","MYC","1","2")

plot <- avr2_yes[,c("dist","MYC_DEX","MYC")]
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
