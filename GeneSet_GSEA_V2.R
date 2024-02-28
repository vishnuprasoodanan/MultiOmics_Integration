install.packages("msigdbr")
library(msigdbr)
library(clusterProfiler)
library(data.table)
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
categories <- msigdbr_collections()
#--------------------------------------------------------------------------------------------------------------
####### function
load_gene_expr <- function(filename){
  genes <- data.frame(fread(filename,sep="\t",head=T), row.names = 1, check.names = F, stringsAsFactors = F)
  genes <- as.matrix(genes)
  
}
#--------------------------------------------------------------------------------------------------------------

## In Rstudio, find the path to the directory where the current script is located.
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
gene_stats <- read.table("results_B_WSDL_HF_Chow.csv", sep = ",", row.names = 1, header = TRUE)
file_names <- list.files(path = current_dir, pattern = "\\.txt$", full.names = TRUE)
#### load data
# Get the names of files ending with .txt in a directory
file_names <- list.files(path = current_dir, pattern = "\\.txt$", full.names = TRUE)

## load gene expression data
for (i in seq_along(file_names)) {
  genes <- load_gene_expr(file_names[i])
  gene_ids_vector <- row.names(genes)
  m_category <- character()
  m_category1 <- character()
  sub_category <- character()
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
  
  # Loop through the rows of the data frame
  for (k in 1:nrow(categories)) {
    # Check if the second column is empty
    if (is.na(categories[k, 2]) || categories[k, 2] == "") {
      # Save the element in the first column to m_category
      m_category <- as.character(categories[k, 1])
      msigdbr_df = msigdbr(species = "mouse", category = m_category)
      msigdbr_t2g = msigdbr_df %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()
      #genes_enrich <- enricher(gene = gene_ids_vector, TERM2GENE = msigdbr_t2g, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
      #ids <- bitr(geneID = gene_ids_vector, fromType = "ENSEMBL", toType  = "UNIPROT", OrgDb = org.Mm.eg.db)
      #geneList_1 = sort(ids$ENTREZID, decreasing = TRUE)
      # y <- gsePathway(geneList_1, 
      #                 pvalueCutoff = 0.2,
      #                 pAdjustMethod = "BH", 
      #                 verbose = FALSE,
      #                 organism = "mouse")
      genes_enrich <- NULL
      
      # Use tryCatch to handle errors
      tryCatch(
        {
          # Your GSEA command
          genes_enrich <- GSEA(ranked_genes, TERM2GENE = msigdbr_t2g)
        },
        error = function(e) {
          # Print the error message or take any other action
          cat("An error occurred:", conditionMessage(e), "\n")
        }
      )
      #genes_enrich <- GSEA(ranked_genes, TERM2GENE = msigdbr_t2g)
      if (!is.null(genes_enrich)) {
        # Write the object to a text file
        output_file2 <- paste0(variables[1], m_category,"_gene_enrich.txt")
        filtered_result <- genes_enrich@result[genes_enrich@result$qvalue < 0.05, ]
        
        # Check if there are any rows left after filtering
        if (nrow(filtered_result) > 0) {
          write.table(filtered_result, file = output_file2, sep = "\t", quote = FALSE)
          print(paste0("Data with qvalue > 0.05 written to the ", output_file2))
          pdf(file = paste0(variables[1], m_category, "_boxplots.pdf"), width = 20, height = 15)
          print(upsetplot(genes_enrich))
          dev.off()
          
        } else {
          print(paste0("No data with qvalue > 0.05 to write.", output_file2))
        }
      }
      else{
        cat("The object is NULL. Moving to the next line of code.\n")
      }
    }
    else {
      # Save elements in first and second columns to m_category1 and sub_category
      m_category1 <- as.character(categories[k, 1])
      sub_category <- as.character(categories[k, 2])
      msigdbr_df = msigdbr(species = "mouse", category = m_category1, subcategory = sub_category)
      msigdbr_t2g = msigdbr_df %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()
      #genes_enrich <- enricher(gene = gene_ids_vector, TERM2GENE = msigdbr_t2g, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
      # Initialize genes_enrich
      genes_enrich <- NULL
      
      # Use tryCatch to handle errors
      tryCatch(
        {
          # Your GSEA command
          genes_enrich <- GSEA(ranked_genes, TERM2GENE = msigdbr_t2g)
        },
        error = function(e) {
          # Print the error message or take any other action
          cat("An error occurred:", conditionMessage(e), "\n")
        }
      )
      #genes_enrich <- GSEA(ranked_genes, TERM2GENE = msigdbr_t2g)
      if (!is.null(genes_enrich)) {
        # Write the object to a text file
        output_file1 <- paste0(variables[1], m_category1, sub_category, "_gene_enrich.txt")
        filtered_result <- genes_enrich@result[genes_enrich@result$qvalue < 0.05, ]
        
        # Check if there are any rows left after filtering
        if (nrow(filtered_result) > 0) {
          write.table(filtered_result, file = output_file1, sep = "\t", quote = FALSE)
          print(paste0("Data with qvalue > 0.05 written to the file.", output_file1))
          pdf(file = paste0(variables[1], m_category1,sub_category, "_boxplots.pdf"), width = 20, height = 15)
          print(upsetplot(genes_enrich))
          dev.off()
        } else {
          print(paste0("No data with qvalue > 0.05 to write.", output_file1))
        }
      }
      else {
        cat("The object is NULL. Moving to the next line of code.\n")
      }
    }
  }
}
#---------------------------------------------------------------------------------------------------------
#Use the gene sets data frame for clusterProfiler with genes as gene symbols.
msigdbr_t2g = msigdbr_df %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
enricher(gene = gene_symbols_vector, TERM2GENE = msigdbr_t2g, ...)

#Use the gene sets data frame for fgsea.
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
fgsea(pathways = msigdbr_list, ...)

#Use the gene sets data frame for GSVA.
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
gsva(gset.idx.list = msigdbr_list, ...)

#----------------------------------------------------------------------------------------------
