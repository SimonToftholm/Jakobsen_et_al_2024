##############################################
################## FIGURE 2 ##################
##############################################

################# Figure 2A ################## 
#The stronger the enhancer are the more regulated while promoters seems less affected by MYC inhibition BT549
library(tidyverse)
BT549_count <- read.delim("MYC_BT549_Consensus.count.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Gene.Name","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_BT549..Tag.Count.in.1000.bp..40939315.0.Total..normalization.factor...0.24..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Gene.Name","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_BT549..Tag.Count.in.1000.bp..40939315.0.Total..normalization.factor...0.24..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh

plot2 <- read.table("DE_analysis_KO_BT549.txt")
plot2 <- plot2[plot2$padj_MYC_vs_Crtl < .05,]
plot2 <- plot2[,c("symbol", "Log2FC_MYC_vs_Crtl")]
plot3 <- left_join(plot, plot2)
plot3[is.na(plot3)] = 0
plot3 <- plot3[!duplicated(plot3$symbol),]

plot4 <- plot3[1:16350,2:5]
plot4 <- plot4[order(plot4$count),]
plot5 <- rowsum(plot4,rep(1:327,each=50))

ggplot() +
  geom_smooth(data= plot5 ,mapping=aes(x=rank(-count),y=Count_enh/50)) +
  geom_smooth(data= plot5 ,mapping=aes(x=rank(-count),y=Count_prom/50)) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + xlab("Ranked Promoter - Enhancer dominace") + ylab("MYC tag density") +
  theme(plot.title = element_text(hjust = 0.5)) 

### Add enhancer connections
gene_enhancer <- read.delim("MYCall_enhancer_gene.txt", header = F)
gene_enhancer <- gene_enhancer[,c("V9","V1")]
gene_enhancer$V1 <- rep(1)
colnames(gene_enhancer) <- c("symbol","test")
test <- left_join(plot3, gene_enhancer)
test[is.na(test)] <- 0
test$scale <- rep("Connection")
test <- test[!duplicated(test$symbol),]
test <- test[test$Log2FC_MYC_vs_Crtl < 0 | test$Log2FC_MYC_vs_Crtl > 0, ]

ggplot(data= test,mapping=aes(x=rank(-count),y=test)) + 
  geom_smooth() + theme(axis.title = element_text(size = 20),
                        axis.text = element_text(size = 20, colour = "black"),
                        panel.background = element_blank(),
                        panel.grid = element_blank(),
                        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + xlab("Ranked Promoter - Enhancer dominace") + ylab("Density of Enh-Pro connections") + coord_cartesian(ylim=c(0,0.3))

################# Figure 2B ################## 
# Gene expression over the time frame
library(mltools)

Up <- plot3[plot3$Log2FC_MYC_vs_Crtl > 0, ]

Up[, "group"] <- bin_data(Up$count, bins=4, binType = "quantile")
one.way <- aov(abs(Log2FC_MYC_vs_Crtl) ~ group, data = Up)
summary(one.way)

ggplot() +
  geom_boxplot(data=Up, mapping=aes(x=group, y=(Log2FC_MYC_vs_Crtl)), alpha=0) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + xlab("Ranked Promoter - Enhancer dominace") + ylab("Gene Log2fc(siMYC/siCrlt)") +
  theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(ylim=c(0,3))

Down <- plot3[plot3$Log2FC_MYC_vs_Crtl < 0 , ]

Down[, "group"] <- bin_data(Down$count, bins=10, binType = "quantile")
one.way <- aov(abs(Log2FC_MYC_vs_Crtl) ~ group, data = Down)
summary(one.way)

ggplot() +
  geom_boxplot(data=Down, mapping=aes(x=group, y=(Log2FC_MYC_vs_Crtl)), alpha=0) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + xlab("Ranked Promoter - Enhancer dominace") + ylab("Gene Log2fc(siMYC/siCrlt)") +
  theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(ylim=c(-1.5,0))

################# Figure 2C ################## 
library(fgsea)
plot3_DEG <- plot3[plot3$Log2FC_MYC_vs_Crtl > 0 | plot3$Log2FC_MYC_vs_Crtl < 0,]

gseaDat <- plot3_DEG
ranks <- gseaDat$count
names(ranks) <- gseaDat$symbol
head(ranks)
barplot(sort(ranks, decreasing = T))

pathways <- gmtPathways("msigdb.v7.4.symbols.gmt")

fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize = 10000, nperm=1000)

head(fgseaRes[order(-NES, abs(NES)), ], n=20)

df <- as.data.frame(fgseaRes)
df_sig <- df[df$pval   < 0.05,]
rownames(df_sig) <- df_sig$pathway
df_plot <- df_sig[grep("HALLMARK_*", df_sig$pathway),]


ggplot(data=df_plot, aes(x=reorder(pathway, NES), y=NES, fill = NES)) +geom_bar(stat="identity", position = "dodge") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom") + scale_fill_viridis_c(option="E") + coord_flip() + xlab(" ") + ylab(" ")

################# Figure 2D ################## 
# MCF7
BT549_count <- read.delim("MYC_MCF7_Consensus.count.txt")
BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_MC7..Tag.Count.in.1000.bp..42607752.5.Total..normalization.factor...0.23..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_MC7..Tag.Count.in.1000.bp..42607752.5.Total..normalization.factor...0.23..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

MCF7_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
MCF7_prom <- as.data.frame(MCF7_prom@result)

MCF7_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
MCF7_enh <- as.data.frame(MCF7_enh@result)



# BT549
BT549_count <- read.delim("MYC_MM1S_enhancer.txt")
BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","MYC_MM1S.Tag.Count.in.1000.bp..19719018.0.Total..normalization.factor...0.51..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","MYC_MM1S.Tag.Count.in.1000.bp..19719018.0.Total..normalization.factor...0.51..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

BT549_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
BT549_prom <- as.data.frame(BT549_prom@result)

BT549_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
BT549_enh <- as.data.frame(BT549_enh@result)




# H2171
BT549_count <- read.delim("MYC_H2171_enhancer.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","MYC_H2171..Tag.Count.in.1000.bp..30710431.0.Total..normalization.factor...0.33..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","MYC_H2171..Tag.Count.in.1000.bp..30710431.0.Total..normalization.factor...0.33..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

H2171_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
H2171_prom <- as.data.frame(H2171_prom@result)

H2171_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
H2171_enh <- as.data.frame(H2171_enh@result)





# U2OS
BT549_count <- read.delim("MYC_U2OS_enhancer.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","MYC_U2OS.Tag.Count.in.1000.bp..5881570.0.Total..normalization.factor...1.70..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","MYC_U2OS.Tag.Count.in.1000.bp..5881570.0.Total..normalization.factor...1.70..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

U2OS_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
U2OS_prom <- as.data.frame(U2OS_prom@result)

U2OS_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
U2OS_enh <- as.data.frame(U2OS_enh@result)


# H128
BT549_count <- read.delim("MYC_H128_enhancer.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","MYC_H128.Tag.Count.in.1000.bp..16989736.0.Total..normalization.factor...0.59..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","MYC_H128.Tag.Count.in.1000.bp..16989736.0.Total..normalization.factor...0.59..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

H128_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
H128_prom <- as.data.frame(H128_prom@result)

H128_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
H128_enh <- as.data.frame(H128_enh@result)

# HT1080
BT549_count <- read.delim("MYC_HT1080_enhancer.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","MYC_HT1080.Tag.Count.in.1000.bp..15495594.0.Total..normalization.factor...0.65..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","MYC_HT1080.Tag.Count.in.1000.bp..15495594.0.Total..normalization.factor...0.65..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

HT1080_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
HT1080_prom <- as.data.frame(HT1080_prom@result)

HT1080_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
HT1080_enh <- as.data.frame(HT1080_enh@result)

# Hela
BT549_count <- read.delim("MYC_HeLA_enhancer.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","MYC_HeLA.Tag.Count.in.1000.bp..12570053.0.Total..normalization.factor...0.80..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","MYC_HeLA.Tag.Count.in.1000.bp..12570053.0.Total..normalization.factor...0.80..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

HeLa_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
HeLa_prom <- as.data.frame(HeLa_prom@result)

HeLa_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
HeLa_enh <- as.data.frame(HeLa_enh@result)

# MM1S
BT549_count <- read.delim("MYC_MM1S_enhancer.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","MYC_MM1S.Tag.Count.in.1000.bp..19719018.0.Total..normalization.factor...0.51..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","MYC_MM1S.Tag.Count.in.1000.bp..19719018.0.Total..normalization.factor...0.51..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

MM1S_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
MM1S_prom <- as.data.frame(MM1S_prom@result)

MM1S_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
MM1S_enh <- as.data.frame(MM1S_enh@result)

# A493
BT549_count <- read.delim("MYC_A493_enhancer.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Entrez.ID","MYC_A493.Tag.Count.in.1000.bp..16697286.0.Total..normalization.factor...0.60..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Entrez.ID","MYC_A493.Tag.Count.in.1000.bp..16697286.0.Total..normalization.factor...0.60..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Entrez.ID) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh
plot <- plot[complete.cases(plot$symbol),]

plot_prom <- plot[plot$count>quantile(plot$count, probs=0.75),"symbol"]
plot_enh <- plot[-plot$count>quantile(-plot$count, probs=0.75),"symbol"]

A493_prom <- enricher(plot_prom$symbol, TERM2GENE=m_t2g)
A493_prom <- as.data.frame(A493_prom@result)

A493_enh <- enricher(plot_enh$symbol, TERM2GENE=m_t2g)
A493_enh <- as.data.frame(A493_enh@result)


# Heatmap
H128_enh <- H128_enh[,c("ID","p.adjust")]
colnames(H128_enh) <- c("pathway","H128")

HT1080_enh <- HT1080_enh[,c("ID","p.adjust")]
colnames(HT1080_enh) <- c("pathway","HT1080")

HeLa_enh <- HeLa_enh[,c("ID","p.adjust")]
colnames(HeLa_enh) <- c("pathway","HeLa")

MM1S_enh <- MM1S_enh[,c("ID","p.adjust")]
colnames(MM1S_enh) <- c("pathway","MM1S")

A493_enh <- A493_enh[,c("ID","p.adjust")]
colnames(A493_enh) <- c("pathway","A493")

BT549_enh <- BT549_enh[,c("ID","p.adjust")]
colnames(BT549_enh) <- c("pathway","BT549")

MCF7_enh <- MCF7_enh[,c("ID","p.adjust")]
colnames(MCF7_enh) <- c("pathway","MCF7")

H2171_enh <- H2171_enh[,c("ID","p.adjust")]
colnames(H2171_enh) <- c("pathway","H2171")

U2OS_enh <- U2OS_enh[,c("ID","p.adjust")]
colnames(U2OS_enh) <- c("pathway","U2OS")


heatmap <- full_join(BT549_enh, MCF7_enh)
heatmap <- full_join(heatmap, H2171_enh)
heatmap <- full_join(heatmap, U2OS_enh )
heatmap <- full_join(heatmap, HT1080_enh)
heatmap <- full_join(heatmap, HeLa_enh)
heatmap <- full_join(heatmap, MM1S_enh)
heatmap <- full_join(heatmap, A493_enh)

rownames(heatmap) <- heatmap$pathway
heatmap <- heatmap[,-1]


heatmap <- heatmap[heatmap$BT549 < 0.05 |
                     heatmap$MCF7 < 0.05 |
                     heatmap$H2171 < 0.05 |
                     heatmap$U2OS < 0.05 |
                     heatmap$HT1080 < 0.05 |
                     heatmap$HeLa < 0.05 |
                     heatmap$MM1S < 0.05 |
                     heatmap$A493 < 0.05 ,]
heatmap <- heatmap[complete.cases(heatmap),]
heatmap[heatmap > 0.05] <- 0.05 
pheatmap(heatmap, scale = "none", color = colorRampPalette(c("#253884","white"))(200))


H128_prom <- H128_prom[,c("ID","p.adjust")]
colnames(H128_prom) <- c("pathway","H128")

HT1080_prom <- HT1080_prom[,c("ID","p.adjust")]
colnames(HT1080_prom) <- c("pathway","HT1080")

HeLa_prom <- HeLa_prom[,c("ID","p.adjust")]
colnames(HeLa_prom) <- c("pathway","HeLa")

MM1S_prom <- MM1S_prom[,c("ID","p.adjust")]
colnames(MM1S_prom) <- c("pathway","MM1S")

A493_prom <- A493_prom[,c("ID","p.adjust")]
colnames(A493_prom) <- c("pathway","A493")

BT549_prom <- BT549_prom[,c("ID","p.adjust")]
colnames(BT549_prom) <- c("pathway","BT549")

MCF7_prom <- MCF7_prom[,c("ID","p.adjust")]
colnames(MCF7_prom) <- c("pathway","MCF7")

H2171_prom <- H2171_prom[,c("ID","p.adjust")]
colnames(H2171_prom) <- c("pathway","H2171")

U2OS_prom <- U2OS_prom[,c("ID","p.adjust")]
colnames(U2OS_prom) <- c("pathway","U2OS")


heatmap <- full_join(BT549_prom, MCF7_prom)
heatmap <- full_join(heatmap, H2171_prom)
heatmap <- full_join(heatmap, U2OS_prom )
heatmap <- full_join(heatmap, HT1080_prom)
heatmap <- full_join(heatmap, HeLa_prom)
heatmap <- full_join(heatmap, MM1S_prom)
heatmap <- full_join(heatmap, A493_prom)


rownames(heatmap) <- heatmap$pathway
heatmap <- heatmap[,-1]


heatmap <- heatmap[heatmap$BT549 < 0.05 |
                     heatmap$MCF7 < 0.05 |
                     heatmap$H2171 < 0.05 |
                     heatmap$U2OS < 0.05 |
                     heatmap$HT1080 < 0.05 |
                     heatmap$HeLa < 0.05 |
                     heatmap$MM1S < 0.05 |
                     heatmap$A493 < 0.05 ,]

heatmap <- heatmap[complete.cases(heatmap),]
heatmap[heatmap > 0.05] <- 0.05
pheatmap(heatmap, scale = "none", color = colorRampPalette(c("#B11518","white"))(200))





################# Figure 2E-F ################## 
library(cBioPortalData)
library(ggplot2)
library(survminer)
library(survival)
library(ggsurvfit)

## Run this if download does not work
mb_exp2 <- read.delim("metabrix_expression.txt", sep = " ")

#Import metastasis data from Rueda et al 2019
relapse<-read.delim("NIHMS1520488-supplement-Supp_Table_6.txt",h=T) #This is clinical data for the metabric cohort downloaded as suppl table S6 from Rueda et al, Nature 2019. The death column is the same as from cbioportal as expected. Info on each column can be found in 41586_2019_1007_MOESM10_ESM_Metabric_variableDescription.xls.

# Divide genes into quantiles
plot4 <- plot3[plot3$Log2FC_MYC_vs_Crtl < -0.5,]
quantile(plot4$count, probs = c(0.3,0.7))

quant_1 <- plot4[plot4$count > 86.100,]
quant_4 <- plot4[plot4$count < -33.782,]

# ER Versus ER+ or HER2+
quant_4_1 <- quant_4
gene_info<-as.data.frame(t(mb_exp2[rownames(mb_exp2) %in% quant_4_1$symbol,])) #Get gene expression for selected genes
gene_info$METABRIC.ID<-rownames(gene_info)
relapse2<-merge(relapse,gene_info,by="METABRIC.ID")

relapse2$SUM_sum <- colMeans((t(relapse2[,36:252]))) 

# High and Low in TNBC
km<-relapse2 %>% filter(Histological.Type=="IDC", ER.Expr=="-",PR.Expr=="-",Her2.Expr=="-")

km$type<-"expression"
km[km$SUM_sum>quantile(km$SUM_sum, probs=0.5),"type"]<-"High" 
km[km$SUM_sum<quantile(km$SUM_sum, probs=0.5),"type"]<-"Low" 

km<-km[km$type!="expression" ,]

km_fit <- survfit(Surv(T, DeathBreast) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = F, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)

km_fit <- survfit(Surv(TDR, DR) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)



# ER Versus ER+ or HER2+
quant_4_1 <- quant_1
gene_info<-as.data.frame(t(mb_exp2[rownames(mb_exp2) %in% quant_4_1$symbol,])) #Get gene expression for selected genes
gene_info$METABRIC.ID<-rownames(gene_info)
relapse2<-merge(relapse,gene_info,by="METABRIC.ID")

relapse2$SUM_sum <- colMeans((t(relapse2[,36:253]))) 

# High and Low in TNBC
km<-relapse2 %>% filter(Histological.Type=="IDC", ER.Expr=="-",PR.Expr=="-",Her2.Expr=="-")

km$type<-"expression"
km[km$SUM_sum>quantile(km$SUM_sum, probs=0.5),"type"]<-"High" 
km[km$SUM_sum<quantile(km$SUM_sum, probs=0.5),"type"]<-"Low" 

km<-km[km$type!="expression" ,]

km_fit <- survfit(Surv(T, DeathBreast) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)


km_fit <- survfit(Surv(TDR, DR) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)


# ER Versus ER+ or HER2+
quant_4_1 <- quant_4
gene_info<-as.data.frame(t(mb_exp2[rownames(mb_exp2) %in% quant_4_1$symbol,])) #Get gene expression for selected genes
gene_info$METABRIC.ID<-rownames(gene_info)
relapse2<-merge(relapse,gene_info,by="METABRIC.ID")

relapse2$SUM_sum <- colMeans((t(relapse2[,36:252]))) 

# High and Low in TNBC
km<-relapse2 %>% filter(ER.Expr=="+", Her2.Expr=="-")
km <- km[complete.cases(km$SUM_sum),]

km$type<-"expression"
km[km$SUM_sum>quantile(km$SUM_sum, probs=0.5),"type"]<-"High" 
km[km$SUM_sum<quantile(km$SUM_sum, probs=0.5),"type"]<-"Low" 

km<-km[km$type!="expression" ,]

km_fit <- survfit(Surv(T, DeathBreast) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)


km_fit <- survfit(Surv(TDR, DR) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)



# ER Versus ER+ or HER2+
quant_4_1 <- quant_1
gene_info<-as.data.frame(t(mb_exp2[rownames(mb_exp2) %in% quant_4_1$symbol,])) #Get gene expression for selected genes
gene_info$METABRIC.ID<-rownames(gene_info)
relapse2<-merge(relapse,gene_info,by="METABRIC.ID")

relapse2$SUM_sum <- colMeans((t(relapse2[,36:253]))) 

# High and Low in TNBC
km<-relapse2 %>% filter(Histological.Type=="IDC", (Pam50Subtype=="LumA"))
km <- km[complete.cases(km$SUM_sum),]

km$type<-"expression"
km[km$SUM_sum>quantile(km$SUM_sum, probs=0.5),"type"]<-"High" 
km[km$SUM_sum<quantile(km$SUM_sum, probs=0.5),"type"]<-"Low" 

km<-km[km$type!="expression" ,]

km_fit <- survfit(Surv(T, DeathBreast) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)


km_fit <- survfit(Surv(TDR, DR) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)


#### MCF7 ####
library(tidyverse)
BT549_count <- read.delim("MYC_MCF7_Consensus.count.txt")

BT549_count$Ann <- sub(" .*","",BT549_count$Annotation)
BT549_anno_Enhancer <- BT549_count[BT549_count$Ann == "Intergenic" | BT549_count$Ann == "intron",]
BT549_anno_Enhancer <- BT549_anno_Enhancer[(BT549_anno_Enhancer$Distance.to.TSS < 50000 & BT549_anno_Enhancer$Detailed.Annotation > 3000) & BT549_anno_Enhancer$Distance.to.TSS > -50000,]
BT549_anno_Promoter <- BT549_count[BT549_count$Ann == "promoter-TSS",]
BT549_anno_Promoter <- BT549_anno_Promoter[!duplicated(BT549_anno_Promoter$Gene.Name),]


BT549_anno_Enhancer_strengh <- BT549_anno_Enhancer[,c("Gene.Name","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_MC7..Tag.Count.in.1000.bp..42607752.5.Total..normalization.factor...0.23..effective.total...10000000.")]
BT549_anno_Enhancer_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> BT549_anno_Enhancer_strengh
colnames(BT549_anno_Enhancer_strengh) <- c("symbol","Count_enh")

BT549_anno_promoter_strengh <- BT549_anno_Promoter[,c("Gene.Name","...MYC_Project.MYC_MED1_ChIP.tagdir.Combined_MYC_MC7..Tag.Count.in.1000.bp..42607752.5.Total..normalization.factor...0.23..effective.total...10000000.")]
BT549_anno_promoter_strengh %>% group_by(Gene.Name) %>% summarise_all(funs(sum)) -> BT549_anno_promoter_strengh
colnames(BT549_anno_promoter_strengh) <- c("symbol","Count_prom")

plot <- full_join(BT549_anno_Enhancer_strengh, BT549_anno_promoter_strengh)
plot[is.na(plot)] = 0
plot$count <- plot$Count_prom-plot$Count_enh

plot3 <- plot[!duplicated(plot$symbol),]

plot4 <- plot3[1:16350,2:4]
plot4 <- plot4[order(plot4$count),]
plot5 <- rowsum(plot4,rep(1:327,each=50))

ggplot() +
  geom_smooth(data= plot5 ,mapping=aes(x=rank(-count),y=Count_enh/50)) +
  geom_smooth(data= plot5 ,mapping=aes(x=rank(-count),y=Count_prom/50)) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + xlab("Ranked Promoter - Enhancer dominace") + ylab("MYC tag density") +
  theme(plot.title = element_text(hjust = 0.5)) 



library(cBioPortalData)
library(ggplot2)
library(survminer)
library(survival)
library(ggsurvfit)
biobrca_mb<-cBioDataPack("brca_metabric", ask = FALSE) #Import everything for metabric
mb_exp<-brca_mb[["expression_median"]]  #Extract mRNA data. The format is: Expression log intensity levels (Illumina Human v3 microarray).
mb_exp2<-as.data.frame(assay(mb_exp)) #Get expression matrix

## Run this if download does not work
mb_exp2 <- read.delim("metabrix_expression.txt", sep = " ")

#Import metastasis data from Rueda et al 2019
relapse<-read.delim("NIHMS1520488-supplement-Supp_Table_6.txt",h=T) #This is clinical data for the metabric cohort downloaded as suppl table S6 from Rueda et al, Nature 2019. The death column is the same as from cbioportal as expected. Info on each column can be found in 41586_2019_1007_MOESM10_ESM_Metabric_variableDescription.xls.

# Divide genes into quantiles
quantile(plot3$count, probs = c(0.1,0.9))

quant_1 <- plot3[plot3$count > 78.4,]
quant_4 <- plot3[plot3$count < -138.4,]

# ER Versus ER+ or HER2+
quant_4_1 <- quant_4
gene_info<-as.data.frame(t(mb_exp2[rownames(mb_exp2) %in% quant_4_1$symbol,])) #Get gene expression for selected genes
gene_info$METABRIC.ID<-rownames(gene_info)
relapse2<-merge(relapse,gene_info,by="METABRIC.ID")

relapse2$SUM_sum <- colMeans(scale(t(relapse2[,36:1104]))) 

# High and Low in TNBC
km<-relapse2 %>% filter(Histological.Type=="IDC", ER.Expr=="-",PR.Expr=="-",Her2.Expr=="-")

km$type<-"expression"
km[km$SUM_sum>quantile(km$SUM_sum, probs=0.8),"type"]<-"High" 
km[km$SUM_sum<quantile(km$SUM_sum, probs=0.2),"type"]<-"Low" 

km<-km[km$type!="expression" ,]

km_fit <- survfit(Surv(T, DeathBreast) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)

km_fit <- survfit(Surv(TDR, DR) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)




# ER Versus ER+ or HER2+
quant_4_1 <- quant_1
gene_info<-as.data.frame(t(mb_exp2[rownames(mb_exp2) %in% quant_4_1$symbol,])) #Get gene expression for selected genes
gene_info$METABRIC.ID<-rownames(gene_info)
relapse2<-merge(relapse,gene_info,by="METABRIC.ID")

relapse2$SUM_sum <- colMeans(scale(t(relapse2[,36:1270]))) 

# High and Low in TNBC
km<-relapse2 %>% filter(Histological.Type=="IDC", ER.Expr=="-",PR.Expr=="-",Her2.Expr=="-")

km$type<-"expression"
km[km$SUM_sum>quantile(km$SUM_sum, probs=0.8),"type"]<-"High" 
km[km$SUM_sum<quantile(km$SUM_sum, probs=0.2),"type"]<-"Low" 

km<-km[km$type!="expression" ,]

km_fit <- survfit(Surv(T, DeathBreast) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)


km_fit <- survfit(Surv(TDR, DR) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)


# ER Versus ER+ or HER2+
quant_4_1 <- quant_4
gene_info<-as.data.frame(t(mb_exp2[rownames(mb_exp2) %in% quant_4_1$symbol,])) #Get gene expression for selected genes
gene_info$METABRIC.ID<-rownames(gene_info)
relapse2<-merge(relapse,gene_info,by="METABRIC.ID")

relapse2$SUM_sum <- colMeans(scale(t(relapse2[,36:1104]))) 

# High and Low in TNBC
km<-relapse2 %>% filter(ER.Expr=="+", Her2.Expr=="-")
km <- km[complete.cases(km$SUM_sum),]

km$type<-"expression"
km[km$SUM_sum>quantile(km$SUM_sum, probs=0.8),"type"]<-"High" 
km[km$SUM_sum<quantile(km$SUM_sum, probs=0.2),"type"]<-"Low" 

km<-km[km$type!="expression" ,]

km_fit <- survfit(Surv(T, DeathBreast) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)


km_fit <- survfit(Surv(TDR, DR) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)



# ER Versus ER+ or HER2+
quant_4_1 <- quant_1
gene_info<-as.data.frame(t(mb_exp2[rownames(mb_exp2) %in% quant_4_1$symbol,])) #Get gene expression for selected genes
gene_info$METABRIC.ID<-rownames(gene_info)
relapse2<-merge(relapse,gene_info,by="METABRIC.ID")

relapse2$SUM_sum <- colMeans(scale(t(relapse2[,36:1270]))) 

# High and Low in TNBC
km<-relapse2 %>% filter(ER.Expr=="+", Her2.Expr=="-")
km <- km[complete.cases(km$SUM_sum),]

km$type<-"expression"
km[km$SUM_sum>quantile(km$SUM_sum, probs=0.8),"type"]<-"High" 
km[km$SUM_sum<quantile(km$SUM_sum, probs=0.2),"type"]<-"Low" 

km<-km[km$type!="expression" ,]

km_fit <- survfit(Surv(T, DeathBreast) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)


km_fit <- survfit(Surv(TDR, DR) ~ type, data=km) #Prepare data for KM plot. Analyse distant metastasis-free survival.
surv_pvalue(km_fit)
ggsurvplot(km_fit,pval = T, pval.method = T, conf.int = F,
           ggtheme = theme(axis.title = element_text(size = 20),
                           axis.text = element_text(size = 20, colour = "black"),
                           panel.background = element_blank(),
                           panel.grid = element_blank(),
                           axis.line = element_line(color = "black"),panel.border = element_rect(colour = "black", fill = NA), legend.position = "bottom"), palette = c("#3F7A13","#828480"), size = 2)



