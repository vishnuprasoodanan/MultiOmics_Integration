if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IHW")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("variancePartition")

library(DESeq2)
library(IHW)
library(ggplot2)

sampletable <- read.table("sample_sheet.txt", header=T, sep="\t")
rownames(sampletable) <- sampletable$SampleName
count_matrix <- as.matrix(read.delim("Kallisto_est_count.tsv", header=T, sep="\t", row.names=1))
count_matrix<-round(count_matrix)
se_star_matrix <- DESeqDataSetFromMatrix(countData = count_matrix,
                                         colData = sampletable,
                                         design = ~gender + Group)
#tx2gene <- read.table("tx2gene.gencode.v29.csv", sep="\t", header=F)

keep <- rowSums(counts(se_star_matrix) >= 30) >= 4
dds <- se_star_matrix[keep,]
nrow(dds)

####variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
star_matrix_vsd_df <- data.frame(assay(vsd))
write.table(star_matrix_vsd_df, file="star_count_matrix_vsd.txt", sep = "\t")

####regularized-logarithm transformation or rlog
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
star_matrix_rld_df <- data.frame(assay(rld))
#####Plot Normalization
library("hexbin")
library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

####### Plot heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDists

library("pheatmap")
library("RColorBrewer")

jpeg("Heatmap_samples.jpg", height = 7, width = 7, units = 'in', res = 600)
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$SampleName, vsd$Group, vsd$MucusGrowth, sep = " - " )
colnames(sampleDistMatrix) <- paste( vsd$SampleName, vsd$Group, vsd$MucusGrowth, sep = " - " )
colors <- colorRampPalette( rev(brewer.pal(9, "BrBG")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()
####### Plot heatmap using Poisson Distance 
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))

jpeg("Heatmap_Poisson.jpg", height = 7, width = 7, units = 'in', res = 600)
samplePoisDistMatrix <- as.matrix( poisd$dd)
rownames(samplePoisDistMatrix) <- paste( vsd$SampleName, vsd$MucusGrowth, vsd$DonorDiet, sep=" - " )
colnames(samplePoisDistMatrix) <- paste( vsd$SampleName, vsd$MucusGrowth, vsd$DonorDiet, sep=" - " )
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)
dev.off()
###### PCA using vst-data
plotPCA(vsd, intgroup = c("MucusGrowth", "DonorDiet"))

jpeg("PCA_vst2.jpg", height = 7, width = 7, units = 'in', res = 600)

pcaData <- plotPCA(vsd, intgroup = c("Group", "MucusGrowth", "MouseDiet"), returnData = TRUE)
pcaData

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = Group, shape = MucusGrowth, label = name, group= MouseDiet)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + theme_bw() + stat_ellipse(geom = "polygon", aes(fill = MouseDiet), alpha = 0.25) +
  ggtitle("PCA with VST data") + geom_text()

dev.off()
###### PCA using rld-data
jpeg("PCA_rld.jpg", height = 7, width = 7, units = 'in', res = 600)

pcaData <- plotPCA(rld, intgroup = c("Group", "MucusGrowth", "DonorDiet"), returnData = TRUE)
pcaData

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = Group, shape = MucusGrowth, label = name, group= DonorDiet)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + theme_bw() + stat_ellipse(geom = "polygon", aes(fill = DonorDiet), alpha = 0.25) +
  ggtitle("PCA with RLD data") + geom_text()

dev.off()

#### generalized principal component analysis, or GLM-PCA

library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$Group <- dds$Group
gpca.dat$MucusGrowth <- dds$MucusGrowth
gpca.dat$DonorDiet <- dds$DonorDiet
gpca.dat$name <- dds$SampleName
gpca.dat$MouseDiet <- dds$MouseDiet
jpeg("generalized_PCA_Group.jpg", height = 7, width = 7, units = 'in', res = 600)
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = Group, shape = MucusGrowth, label = name, group= Group)) + 
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA") + theme_bw() +geom_text() + stat_ellipse(geom = "polygon", aes(fill = Group), alpha = 0.25)
dev.off()
####################################
#---- Multidimensional scaling (MDS)
###################################
#on vsd data
sampleDists <- dist(t(assay(vsd)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
jpeg("MSD_plot_vsd_MucusGrowth.jpg", height = 7, width = 7, units = 'in', res = 600)
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = Group, shape = MucusGrowth, label = SampleName, group= MucusGrowth)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data") + theme_bw()+ geom_text() + stat_ellipse(geom = "polygon", aes(fill = MucusGrowth), alpha = 0.25)
dev.off()

sampleDists <- dist(t(assay(dds)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
jpeg("MSD_plot_dds_MucusGrowth.jpg", height = 7, width = 7, units = 'in', res = 600)
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = Group, shape = MucusGrowth, label = SampleName, group= MucusGrowth)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS") + theme_bw()+ geom_text() + stat_ellipse(geom = "polygon", aes(fill = MucusGrowth), alpha = 0.25)
dev.off()

#based on Poisson distance on count data
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd)
jpeg("MSD_plot_PoissonDistance.jpg", height = 7, width = 7, units = 'in', res = 600)
mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = MucusGrowth, shape = DonorDiet, label = SampleName)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances") + theme_bw() + geom_text()
dev.off()

#################################
#---- Variance Partition Analysis
##################################
# on vsd data
Genenames<-row.names(star_matrix_vsd_df)
library(variancePartition)
model <- ~ (1|Gender) + (1|DonorDiet) + (1|MucusGrowth) + (1|Group) + (1|MouseDiet)
varPart_model <- fitExtractVarPartModel(star_matrix_vsd_df, model, sampletable)
dim(varPart_model)
save(varPart_model, file="varPart_model.Rdata")
load ("varPart_model.Rdata")
rownames(varPart_model) = make.names(Genenames, unique=TRUE)
varPart_model_df <- data.frame(varPart_model)
write.table(varPart_model_df, file="varPart_model.txt", sep = "\t")
#sort the terms of the model by average variance explained across all genes, so when we plot they will be sorted by overall importance:
vp2 <- sortCols( varPart_model )
# Violin plot
jpeg("variance_partision_vst.jpg", height = 7, width = 7, units = 'in', res = 600)
plotVarPart( vp2  ,label.angle = 90)
dev.off ()

# on raw count data
Genenames<-row.names(star_matrix_rld_df)
library(variancePartition)
model <- ~ (1|Gender) + (1|DonorDiet) + (1|MucusGrowth) + (1|Group) + (1|MouseDiet)
varPart_model <- fitExtractVarPartModel(star_matrix_rld_df, model, sampletable)
dim(varPart_model)
save(varPart_model, file="varPart_model_rawcount.Rdata")
load ("varPart_model_rawcount.Rdata")
rownames(varPart_model) = make.names(Genenames, unique=TRUE)
varPart_model_df <- data.frame(varPart_model)
write.table(varPart_model_df, file="varPart_model_rld.txt", sep = "\t")
#sort the terms of the model by average variance explained across all genes, so when we plot they will be sorted by overall importance:
vp2 <- sortCols( varPart_model )
# Violin plot
jpeg("variance_partision_on_rld.jpg", height = 7, width = 7, units = 'in', res = 600)
plotVarPart( vp2  ,label.angle = 90)
dev.off ()

#### Differential expression
dds <- DESeq(dds)

res <- results(dds, contrast=c("Group","BL_Chow","HF_Chow"), filterFun=ihw)
mcols(res, use.names = TRUE)
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
write.csv(res, file = "results_BL_HF_Chow.csv")

sum(res$pvalue < 0.05, na.rm=TRUE)
sum(res$pvalue < 0.01, na.rm=TRUE)

sum(!is.na(res$pvalue))

sum(res$padj < 0.1, na.rm=TRUE)


resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])


topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("MucusGrowth"))
###################################
# Plot selected genes
###################################
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("MucusGrowth","DonorDiet"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = MucusGrowth, y = count, color = MucusGrowth)) +
  scale_y_log10() +  geom_beeswarm(cex = 5)+ theme_bw()

######

library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)


mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("MucusGrowth","DonorDiet")])
pheatmap(mat, annotation_col = anno)

######

library("IHW")
res.ihw <- results(dds, filterFun=ihw)

###### Annotation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Mm.eg.db")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")

library("AnnotationDbi")
library("org.Mm.eg.db")

columns(org.Mm.eg.db)
ens.str <- rownames(res)
res$symbol <- mapIds(org.Mm.eg.db,
                     keys=ens.str,
                     column=c('SYMBOL'),
                     keytype='ENSEMBL',
                     multiVals="first")
res$entrez <- mapIds(org.Mm.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

resOrdered <- res[order(res$pvalue),]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)
write.csv(resOrderedDF, file = "results.csv")

####volcano plot
library(EnhancedVolcano)
library(DESeq2)
library(IHW)
library(ggplot2)

sampletable <- read.table("sample_sheet.txt", header=T, sep="\t")
rownames(sampletable) <- sampletable$SampleName
count_matrix <- as.matrix(read.delim("STAR_gene_count.tsv", header=T, sep="\t", row.names=1))
count_matrix<-round(count_matrix)
se_star_matrix <- DESeqDataSetFromMatrix(countData = count_matrix,
                                         colData = sampletable,
                                         design = ~Group)

keep <- rowSums(counts(se_star_matrix) >= 30) >= 4
dds <- se_star_matrix[keep,]
nrow(dds)

vsd <- vst(dds, blind = FALSE)
dds <- DESeq(dds)
res <- results(dds, contrast=c("Group","HF_Chow","HF_WSD"), filterFun=ihw)
jpeg("VolcanoPlot_HF_Chow_WSD.jpg", height = 20, width = 20, units = 'in', res = 600)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                pCutoff = 0.05,
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)
dev.off()

pdf(file="VolcanoPlot_HF_Chow_WSD.pdf", height = 20, width = 20)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                pCutoff = 0.05,
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)
dev.off()


library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)


mat  <- assay(vsd)[ topVarGenes, ]
jpeg("Selected_genes_heatmap.jpg", height = 7, width = 7, units = 'in', res = 600)
mat <- as.matrix(read.table(file = "Selected_genes_count_matrix_vsd_86503.txt", header = TRUE, sep = "\t", row.names = 1))
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("MucusGrowth","DonorDiet")])
pheatmap(mat, annotation_col = anno)
dev.off()


library("ggbeeswarm")
library("gridExtra")
pdf(file="gene2.pdf")
par(mfrow=c(2,2))
Colors <- c("red2", "green3")
Colors1 <- c("red4", "green4")
gene_list <- c("ENSMUSG00000066861", "ENSMUSG00000068452", "ENSMUSG00000034459", "ENSMUSG00000025810", "ENSMUSG00000052776", "ENSMUSG00000035692", "ENSMUSG00000042312", "ENSMUSG00000079339", "ENSMUSG00000020641", "ENSMUSG00000076512", "ENSMUSG00000032661", "ENSMUSG00000033355", "ENSMUSG00000029561", "ENSMUSG00000096215", "ENSMUSG00000074896", "ENSMUSG00000033213", "ENSMUSG00000064339", "ENSMUSG00000037921", "ENSMUSG00000038146", "ENSMUSG00000021208", "ENSMUSG00000096001", "ENSMUSG00000041827", "ENSMUSG00000057346", "ENSMUSG00000025498", "ENSMUSG00000027514", "ENSMUSG00000028037", "ENSMUSG00000041134", "ENSMUSG00000064337", "ENSMUSG00000046814", "ENSMUSG00000014313", "ENSMUSG00000027317", "ENSMUSG00000075602", "ENSMUSG00000022371", "ENSMUSG00000030762", "ENSMUSG00000040483", "ENSMUSG00000024168", "ENSMUSG00000027078", "ENSMUSG00000031104", "ENSMUSG00000100975", "ENSMUSG00000050777", "ENSMUSG00000074445", "ENSMUSG00000027225", "ENSMUSG00000068246", "ENSMUSG00000025161", "ENSMUSG00000086503", "ENSMUSG00000032690")
plot_lst <- vector("list", length = length(gene_list))

for (i in 1:length(gene_list)){
  topGene <- gene_list[i]
  
  geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("MucusGrowth","DonorDiet"),
                           returnData = TRUE)
  data1 <- geneCounts
  #P <- ggplot(geneCounts, aes(x = MucusGrowth, y = count, color = MucusGrowth)) +
  #  scale_y_log10() +  geom_beeswarm(cex = 10)+ theme_bw() + ggtitle(topGene)
  P<- ggplot(data1, aes(MucusGrowth, count, fill=MucusGrowth))+
    ggtitle(topGene)+
    labs(y = "count")+
    geom_boxplot(outlier.shape=NA)+ 
    geom_point(aes(colour = factor(data1$MucusGrowth)), position=position_jitterdodge(jitter.width = 0.5))+scale_y_log10()+
    scale_color_manual(values = Colors1)+
    scale_fill_manual(values = Colors)+ theme_bw() + theme(legend.position="none")
  
  plot_lst[[i]] <- P
}
ml <- marrangeGrob(plot_lst, nrow = 2, ncol = 2)
print(ml)
dev.off()
