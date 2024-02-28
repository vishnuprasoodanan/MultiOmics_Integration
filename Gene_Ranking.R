library(fgsea)
library(data.table)
library(ggplot2)
library("AnnotationDbi")
library("org.Mm.eg.db")
library(dplyr)
library(ReactomePA)
library(gridExtra)
library(ggplot2)
library(ggplotify)
library(enrichplot)
set.seed(42)

#--------------------------------------------------------------------------------------------------------------
####### function
load_gene_expr <- function(filename){
  genes <- data.frame(fread(filename,sep="\t",head=T), row.names = 1, check.names = F, stringsAsFactors = F)
  genes <- as.matrix(genes)
  
}
#----------------------------------------------------------------------
#categories <- msigdbr_collections()
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
gene_stats <- read.table("results_B_WSDL_HF_Chow.csv", sep = ",", row.names = 1, header = TRUE)
file_names <- list.files(path = current_dir, pattern = "\\.txt$", full.names = TRUE)

## load gene expression data
for (i in seq_along(file_names)) {
  genes <- load_gene_expr(file_names[i])
  variables <- strsplit(basename(file_names[i]), '_')[[1]]
  sel_genes <- row.names(genes)
  # Making a subset
  subset <- gene_stats[row.names(gene_stats) %in% sel_genes, ]
  subset_sort <- subset[order(-subset$stat), ]
  ranked_genes <- subset_sort$stat
  names(ranked_genes) <- rownames(subset_sort)
  check1 <- as.data.frame(ranked_genes)
  check2<- row.names(check1)

  #columns(org.Hs.eg.db) # returns list of available keytypes
  check3 <- as.data.frame(mapIds(org.Mm.eg.db, keys=check2, column="ENTREZID", keytype="ENSEMBL", multiVals="first"))
  colnames(check3)[1] <- "ENTREZ_ID"
  merged_df <- merge(check1, check3, by = "row.names", all.x = TRUE)
  merged_df_cleaned <- na.omit(merged_df)
  result_df <- as.data.frame(merged_df_cleaned %>% group_by(ENTREZ_ID) %>% sample_n(1))
  rownames(result_df) <- result_df$ENTREZ_ID
  result_df$Row.names <- NULL
  result_df$ENTREZ_ID <- NULL
  merged_df_cleaned_num <- as.vector(result_df$ranked_genes)
  names(merged_df_cleaned_num) <- rownames(result_df)
  merged_df_cleaned_num <- merged_df_cleaned_num[order(-merged_df_cleaned_num)]
  
  
  # gse <- gseMKEGG(geneList = merged_df_cleaned_num,
  #                     organism = 'mmu',
  #                     pvalueCutoff = 1)
  gse <- gsePathway(merged_df_cleaned_num, organism = 'mouse',
             pvalueCutoff = 0.2,
             pAdjustMethod = "BH", 
             verbose = FALSE, eps = 0)
  if (!is.null(gse)) {
    # Write the object to a text file
    output_file2 <- paste0(variables[1], "KEGG_gene_enrich.txt")
    filtered_result <- gse@result[gse@result$pvalue < 1, ]
    # Check if there are any rows left after filtering
    if (nrow(filtered_result) > 0) {
      write.table(filtered_result, file = output_file2, sep = "\t", quote = FALSE)
      print(paste0("Data with qvalue > 0.05 written to the ", output_file2))
      pdf(file = paste0(variables[1], "_boxplots.pdf"), width = 20, height = 15)
      print(upsetplot(gse))
      dev.off()
      # pdf(file = paste0(variables[1], "_gsea_plot.pdf"), width = 20, height = 15)
      # print(gseaplot2(gse, geneSetID = 1:length(rownames(gse@result)), title = variables[1]))
      # dev.off()
      pdf(file = paste0(variables[1], "_gsea_plot.pdf"), width = 20, height = 15)
      plot_lst <- list("list", length = as.numeric(length(rownames(gse@result))))
      for (i in 1:length(rownames(gse@result))) {
        P <- as.ggplot(gseaplot2(gse, geneSetID = i, title = gse$Description[i]))
        plot_lst[[i]] <- P
      }
      ml <- marrangeGrob(plot_lst, nrow = 1, ncol = 1)
      print(ml)
      dev.off()
      
      } else
        {
          print(paste0("No data with qvalue > 0.05 to write.", output_file2))
        }
    }
  else{
    cat("The object is NULL. Moving to the next line of code.\n")
  }
}
# If you want to remove the redundant row names column
rownames(merged_df_cleaned) <- merged_df_cleaned$ENTREZ_ID
merged_df_cleaned$Row.names <- NULL
merged_df_cleaned$ENTREZ_ID <- NULL
# Convert the single-column dataframe to a named vector
merged_df_cleaned_num <- as.vector(merged_df_cleaned$ranked_genes)
names(merged_df_cleaned_num) <- rownames(merged_df_cleaned)
merged_df_cleaned_num <- merged_df_cleaned_num[order(-merged_df_cleaned_num)]




ids <- bitr(row.names(check1), fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")