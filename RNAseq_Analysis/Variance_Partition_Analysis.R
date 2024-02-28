library(DESeq2)
library(IHW)
library(ggplot2)

sampletable <- read.table("sample_sheet_v6.txt", header=T, sep="\t")
rownames(sampletable) <- sampletable$SampleName
count_matrix <- as.matrix(read.delim("STAR_gene_count.tsv", header=T, sep="\t", row.names=1))
count_matrix<-round(count_matrix)
se_star_matrix <- DESeqDataSetFromMatrix(countData = count_matrix,
                                         colData = sampletable,
                                         design = ~ Gender)
#tx2gene <- read.table("tx2gene.gencode.v29.csv", sep="\t", header=F)

keep <- rowSums(counts(se_star_matrix) >= 30) >= 4
dds <- se_star_matrix[keep,]
nrow(dds)

vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
star_matrix_vsd_df <- data.frame(assay(vsd))
write.table(star_matrix_vsd_df, file="star_count_matrix_vsd.txt", sep = "\t")

#################################
#---- Variance Partition Analysis
##################################
# on vsd data
Genenames<-row.names(star_matrix_vsd_df)
library(variancePartition)
#model <- ~ (1|Gender) + (1|DonorDiet) + (1|MucusGrowth) + (1|Group) + (1|MouseDiet) + (1|B_cell) + (1|Cardiac_Muscle) + (1|Epidermis) + (1|Epithelial_cell) + (1|Erythrocyte) + (1|Granulocyte) + (1|Hepatocyte) + (1|Monocyte) + (1|NK_cell) + (1|Stromal_cell) + (1|T_cell)
#model <- ~ (1|Gender) + (1|DonorDiet) + (1|MucusGrowth) + (1|Group) + (1|MouseDiet) + B_cell + Cardiac_Muscle + Epidermis + Epithelial_cell + Erythrocyte + Granulocyte + Hepatocyte + Monocyte + NK_cell + Stromal_cell + T_cell
#model <- ~ (1|Gender) + (1|DonorDiet) + (1|MucusGrowth) + (1|Group) + (1|MouseDiet) + B_cell + Epithelial_cell + T_cell
model <- ~ (1|Gender) + (1|DonorDiet) + (1|MucusGrowth_discrete) + (1|Group) + (1|MouseDiet) + B_cell + Epithelial_cell + T_cell + intestine_large + intestine_small + T.cells_CD8. + T.cells_CD4. + MucusGrowth_continuous
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
jpeg("variance_partision_vst_v4.jpg", height = 7, width = 7, units = 'in', res = 600)
plotVarPart( vp2  ,label.angle = 90)
dev.off ()

######################################
library("genefilter")
topVarGenes <- c("ENSMUSG00000056673", "ENSMUSG00000086503", "ENSMUSG00000068457", "ENSMUSG00000069049", "ENSMUSG00000069045", "ENSMUSG00000099876", "ENSMUSG00000101036", "ENSMUSG00000037369", "ENSMUSG00000025332", "ENSMUSG00000067149", "ENSMUSG00000095079", "ENSMUSG00000041420", "ENSMUSG00000076609", "ENSMUSG00000031226", "ENSMUSG00000035150", "ENSMUSG00000087174", "ENSMUSG00000025209", "ENSMUSG00000030823", "ENSMUSG00000079105", "ENSMUSG00000037033", "ENSMUSG00000039114", "ENSMUSG00000044702", "ENSMUSG00000049112", "ENSMUSG00000029710", "ENSMUSG00000025507", "ENSMUSG00000036529", "ENSMUSG00000076545", "ENSMUSG00000009092", "ENSMUSG00000047632", "ENSMUSG00000047832", "ENSMUSG00000032555", "ENSMUSG00000021175", "ENSMUSG00000042195", "ENSMUSG00000059447", "ENSMUSG00000003435", "ENSMUSG00000028664", "ENSMUSG00000022584", "ENSMUSG00000023110", "ENSMUSG00000094546", "ENSMUSG00000009112", "ENSMUSG00000105703", "ENSMUSG00000024182", "ENSMUSG00000000384", "ENSMUSG00000050097", "ENSMUSG00000022283", "ENSMUSG00000039308", "ENSMUSG00000035206", "ENSMUSG00000097615", "ENSMUSG00000002489", "ENSMUSG00000031171")
mat  <- assay(vsd)[ topVarGenes, ]
jpeg("Selected_genes_heatmap.jpg", height = 7, width = 7, units = 'in', res = 600)
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Gender","DonorDiet")])
pheatmap(mat, annotation_col = anno)
dev.off()

#####################################

