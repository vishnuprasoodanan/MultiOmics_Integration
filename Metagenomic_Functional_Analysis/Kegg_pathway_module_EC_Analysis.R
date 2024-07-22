# Load required libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(labdsv)
library(RColorBrewer)
library(vegan)
library(ape)
library(vegan)
library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(phyloseq)
library(reshape2)
library(ape)
library(mia)
library(Boruta)
library(ANCOMBC)
# Define parameters
Colors <- c("darkolivegreen4", "red", "darkgreen", "salmon4")
alpha_value <- 0.5
circle_alpha <- 0.8

# Function to normalize each row by dividing by the row sum
rel_abundance <- function(df) {
  row_sums <- rowSums(df)
  normalized_df <- sweep(df, 1, row_sums, FUN = "/")
  return(normalized_df)
}

# Function to process the dataframe as specified
process_dataframe <- function(df) {
  # Replace values less than 10 with 0
  df[df < 100] <- 0
  
  # Include rows with fewer than 8 zeros
  row_zero_counts <- rowSums(df == 0)
  filtered_df <- df[row_zero_counts < 8, ]
  
  return(filtered_df)
}

# Function to process the dataframe as specified
process_rpkm_dataframe <- function(df) {
  # Replace values less than 10 with 0
  df[df < 0.5] <- 0
  
  # Include rows with fewer than 8 zeros
  row_zero_counts <- rowSums(df == 0)
  filtered_df <- df[row_zero_counts < 8, ]
  
  return(filtered_df)
}
# Function to normalize each column log2 transformation
log_transform <- function(data) {
  logcounts <- log2(data + 1)
  logcounts_matrix <- as.data.frame(t(logcounts))
  return(logcounts_matrix)
}

# Function to perform CLR transformation
clr_transformation <- function(counts, samples) {
  # Extract row names of df1
  count_row_names <- row.names(counts)
  
  # Create new dataframe with one column 'Taxon' containing row names as values
  tax <- data.frame(Taxon = count_row_names)
  
  # Set row names of df2 to be the same as the values in 'Taxon' column
  row.names(tax) <- tax$Taxon
  
  # Create SummarizedExperiment object
  se <- SummarizedExperiment(assays = list(counts = counts),
                             colData = samples,
                             rowData = tax)
  tse <- as(se, "TreeSummarizedExperiment")
  
  # CLR transformation
  tse <- transformAssay(tse, method = "clr", pseudocount = 1)
  clr_assay <- assays(tse)$clr
  clr_assay_df <- as.data.frame(clr_assay)
  
  return(clr_assay_df)
}

# Function to check if an element is not empty
is_not_empty <- function(x) {
  !is.null(x) && length(x) > 0
}

# Function to move the 'Status' column to the first position in the dataframe
move_status_to_first <- function(df) {
  n_cols <- ncol(df)
  df <- df[, c(n_cols, 1:(n_cols - 1))]
  return(df)
}

# Function to subset a dataframe if columns and rows are provided as vectors
subset_dataframe <- function(df, selected_columns, selected_rows) {
  # Check if selected_columns is empty
  if (length(selected_columns) == 1 && selected_columns == "") {
    selected_columns <- colnames(df)
  } else {
    # Exclude columns that are not in the dataframe
    selected_columns <- selected_columns[selected_columns %in% colnames(df)]
  }
  
  # Check if selected_rows is empty
  if (length(selected_rows) == 1 && selected_rows == "") {
    selected_rows <- 1:nrow(df)
  } else {
    # Exclude rows that are not in the dataframe
    selected_rows <- selected_rows[selected_rows %in% rownames(df)]
  }
  
  # Subset the dataframe
  df_subset <- df[selected_rows, selected_columns, drop = FALSE]
  
  return(df_subset)
}

#Function to merge two dataframes by rownames and add firstcolumn of df2 as first column of merged dataframe
merge_by_rownames <- function(df1, df2) {
  # Ensure rownames are set correctly
  if (is.null(rownames(df1)) || is.null(rownames(df2))) {
    stop("Both dataframes must have rownames.")
  }
  
  # Merge dataframes by rownames
  common_rownames <- intersect(rownames(df1), rownames(df2))
  
  # Subset dataframes by common rownames
  df1_subset <- df1[common_rownames, , drop = FALSE]
  df2_subset <- df2[common_rownames, , drop = FALSE]
  
  # Add the first column of df2 to df1 as the first column
  df1_with_df2_first_col <- cbind(df2_subset[, 1, drop = FALSE], df1_subset)
  
  return(df1_with_df2_first_col)
}

# Split dataframe containing 'Status' in the first column as metadata, taxon, and count dataframes
split_dataframe <- function(df) {
  # Check if dataframe has at least one column
  if (ncol(df) < 1) {
    stop("Dataframe must have at least one column.")
  }
  
  # Part 1: Dataframe containing rownames and the first column
  sub_metadata <- data.frame(Status = df[, 1], row.names = rownames(df))
  
  # Part 2: Dataframe except the first column
  sub_count <- df[, -1, drop = FALSE]
  
  # Part 3: Dataframe with column names as rownames and same column names in the first column as 'Taxon'
  sub_taxa <- data.frame(Taxon = colnames(df), row.names = colnames(df))
  
  # Remove the first row of sub_taxa
  sub_taxa <- sub_taxa[-1, , drop = FALSE]
  
  return(list(sub_metadata = sub_metadata, sub_count = sub_count, sub_taxa = sub_taxa))
}

# Split dataframe containing 'Status' in the first column as metadata, taxon, and count dataframes
split_status_column <- function(df) {
  # Check if the first column is named 'Status'
  if (colnames(df)[1] != "Status") {
    stop("Error: The first column is not named 'Status'")
  }
  
  # Split the values in the 'Status' column based on '_'
  split_values <- strsplit(as.character(df$Status), "_")
  
  # Extract the last value after splitting
  last_values <- sapply(split_values, function(x) tail(x, n=1))
  
  # Extract the second last value after splitting
  second_last_values <- sapply(split_values, function(x) ifelse(length(x) > 1, x[length(x) - 1], NA))
  
  # Create the two new dataframes
  df_last <- df
  df_last$Status <- last_values
  
  df_second_last <- df
  df_second_last$Status <- second_last_values
  
  # Return the two dataframes
  return(list(df_last = df_last, df_second_last = df_second_last))
}

#Function to create boxplots for each dataframe and save as PDF
create_boxplots <- function(df, filename) {
  if (nrow(df) == 0) {
    cat(paste("Dataframe", filename, "has zero elements. Skipping...\n"))
    return(NULL)
  }
  
  plots <- lapply(colnames(df)[-1], function(col) {
    p <- ggplot(df, aes(x = Status, y = !!sym(col), fill = Status)) +
      geom_boxplot(alpha = alpha_value, color = "black") +
      geom_jitter(shape = 16, size = 2, alpha = circle_alpha, aes(color = Status)) +
      scale_color_manual(values = Colors) +
      scale_fill_manual(values = Colors) +
      labs(title = col) +
      theme_minimal() +
      theme(legend.position = "none")
    return(p)
  })
  
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0) {
    cat(paste("No valid plots to create for", filename, ". Skipping...\n"))
    return(NULL)
  }
  
  pdf_file <- file.path(paste0(filename, ".pdf"))
  pdf(pdf_file, width = 12, height = 12)
  
  num_plots <- length(plots)
  num_pages <- ceiling(num_plots / 9)
  
  for (page in seq_len(num_pages)) {
    start_plot <- (page - 1) * 9 + 1
    end_plot <- min(page * 9, num_plots)
    
    if (start_plot <= end_plot) {
      plots_subset <- plots[start_plot:end_plot]
      
      if (length(plots_subset) < 9) {
        num_cols <- 3
        num_rows <- 3
        
        if (length(plots_subset) < 4) {
          num_cols <- 3
          num_rows <- 3
          
          if (length(plots_subset) < 2) {
            num_cols <- 3
            num_rows <- 3
          }
        }
        
        empty_plots <- (num_cols * num_rows) - length(plots_subset)
        for (i in 1:empty_plots) {
          plots_subset <- c(plots_subset, NULL)
        }
      } else {
        num_cols <- 3
        num_rows <- 3
      }
      
      plot_grid <- do.call(grid.arrange, c(plots_subset, ncol = num_cols, nrow = num_rows))
      print(plot_grid)
    }
  }
  
  dev.off()
  cat(paste("PDF file", pdf_file, "created successfully.\n"))
}

# Function to plot PCoA using Bray-Curtis distance and perform PERMANOVA
plot_pcoa_and_permanova <- function(df, output_filename) {
  class <- df$Status
  
  bray_distances <- vegdist(df[, -1], method = "bray")
  bray_pcoa <- pcoa(bray_distances)
  
  mds_var_per <- round(bray_pcoa$values$Relative_eig[1:2] * 100, 1)
  bray_pcoa_matrix <- as.data.frame(bray_pcoa$vectors[, 1:2])
  
  adonis2_result <- adonis2(bray_distances ~ df$Status)
  print(adonis2_result)
  
  p_value <- adonis2_result$`Pr(>F)`[1]
  r_squared <- adonis2_result$R2[1]
  f_value <- adonis2_result$F[1]
  
  bray_pcoa_matrix_new <- cbind(bray_pcoa_matrix, class)
  
  pdf_file <- paste0(output_filename, "PCoA_Plot.pdf")
  pdf(pdf_file, height = 10, width = 10)
  
  plot(bray_pcoa_matrix_new[, 1:2], bg = Colors[as.factor(bray_pcoa_matrix_new$class)], pch = 21, cex = 2,
       xlab = paste0("PCoA", 1, " (", mds_var_per[1], "%)"),
       ylab = paste0("PCoA", 2, " (", mds_var_per[2], "%)"))
  
  ordiellipse(bray_pcoa_matrix_new[, 1:2], bray_pcoa_matrix_new$class, kind = "sd", lwd = 1, lty = 3, draw = "polygon", alpha = 70, col = Colors)
  
  legend("topright", legend = levels(as.factor(df$Status)), col = Colors, lty = 1, cex = 0.7, bg = "white", bty = "n")
  
  text(x = min(bray_pcoa_matrix_new[, 1]), y = max(bray_pcoa_matrix_new[, 2]),
       labels = paste("p-value: ", format(p_value, digits = 3),
                      "\nR-squared: ", format(r_squared, digits = 3),
                      "\nF-value: ", format(f_value, digits = 3)),
       adj = c(0, 1), cex = 0.8, pos = 4)
  
  dev.off()
  
  cat(paste("PCoA plot saved as", pdf_file, "\n"))
}

# Function to create PCoA plot using Euclidean distance and perform PERMANOVA
create_clr_pcoa_plot <- function(clr_assay_df, samples, Colors, output_file) {
  # Calculate Euclidean distances
  euclidean_dist <- vegdist(as.matrix(clr_assay_df), method = "euclidean")
  
  # PCoA analysis
  Uni_pcoa <- pcoa(as.matrix(euclidean_dist))
  mds.var.per <- round(Uni_pcoa$values$Eigenvalues / sum(Uni_pcoa$values$Eigenvalues) * 100, 1)
  pc <- c(1, 2)
  
  # PERMANOVA analysis
  adonis_result <- adonis2(as.matrix(euclidean_dist) ~ samples$Status)
  r_squared <- adonis_result$R2[1]
  p_value <- adonis_result$Pr[1]
  f_value <- adonis_result$F[1]
  
  # Create PCoA plot and save to PDF
  pdf(output_file, height = 10, width = 10)
  plot(Uni_pcoa$vectors[, 1:2], 
       bg = Colors[as.factor(samples$Status)], 
       pch = 21, cex = 2, 
       xlab = paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), 
       ylab = paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
  ordiellipse(Uni_pcoa$vectors[, 1:2], samples$Status, kind = "sd", lwd = 1, lty = 3, 
              draw = "polygon", alpha = 50, col = Colors)
  abline(h = 0, v = 0, col = "gray60")
  
  # Add R-squared, p-value, and F-value to the plot
  legend("topright", legend = paste("R^2 =", round(r_squared, 3), 
                                    "\nP =", round(p_value, 3), 
                                    "\nF =", round(f_value, 3)), 
         bty = "n")
  dev.off()
}

# Function to perform Mantel test and store results
perform_mantel_test <- function(df_list) {
  results <- data.frame(
    dataframe1 = character(),
    dataframe2 = character(),
    Mantel_statistic_r = numeric(),
    Significance = numeric(),
    stringsAsFactors = FALSE
  )
  
  combinations <- combn(names(df_list), 2, simplify = FALSE)
  
  for (pair in combinations) {
    mantel_result <- mantel(df_list[[pair[1]]], df_list[[pair[2]]], method = "spearman", permutations = 999)
    results <- rbind(results, data.frame(
      dataframe1 = pair[1],
      dataframe2 = pair[2],
      Mantel_statistic_r = mantel_result$statistic,
      Significance = mantel_result$signif,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}
#------------------------------------Main code starts here----------------------------------------

# Create output directories if they don't exist
#if (!dir.exists("output_tables")) dir.create("output_tables")
#if (!dir.exists("output_plots")) dir.create("output_plots")

# Read data
#gene_count <- read.csv("Final_GeneCount_100_criteria.txt", sep = "\t", header = TRUE, row.names = 1)
#cazy_count <- read.csv("CAZy_genes_final.txt_summary.tsv_GeneCount.txt", sep = "\t", header = TRUE, row.names = 1)
#ko_count <- read.csv("KEGG_ko_final.txt_summary.tsv_GeneCount.txt", sep = "\t", header = TRUE, row.names = 1)
ec_count <- read.csv("Selected_EC_GeneCount.txt", sep = "\t", header = TRUE, row.names = 1)
ko_module <- read.csv("Selected_KEGG_Module_GeneCount.txt", sep = "\t", header = TRUE, row.names = 1)
ko_path <- read.csv("Selected_Pathways_GeneCount.txt", sep = "\t", header = TRUE, row.names = 1)
#ko_reaction <- read.csv("KEGG_Reaction_final.txt_summary.tsv_GeneCount.txt", sep = "\t", header = TRUE, row.names = 1)
#ko_rclass <- read.csv("KEGG_rclass_genes_final.txt_summary.tsv_GeneCount.txt", sep = "\t", header = TRUE, row.names = 1)

metadata <- read.csv("Metadata.txt", sep = "\t", header = TRUE, row.names = 1)
#cazy_rpkm <- read.csv("HUMAN_CAYMAN_CAZy_RPKM.txt", sep = "\t", header = TRUE, row.names = 1)

fun_count_list <- list(#gene_count = gene_count,
                       #cazy_count = cazy_count,
                       #ko_count = ko_count,
                       ec_count = ec_count,
                       ko_module = ko_module,
                       ko_path = ko_path
                       #ko_reaction = ko_reaction,
                       #ko_rclass = ko_rclass)
                      )
#rm(gene_count, cazy_count, ko_count, ec_count, ko_module, ko_path, ko_reaction, ko_rclass)

# List to store the processed and CLR-transformed dataframes
fun_clr_list <- list()
fun_log2_list <- list()
fun_rel_list <- list()
# Loop through each dataframe in fun_count_list
for (df_name in names(fun_count_list)) {
  df <- fun_count_list[[df_name]]
  
  # Process the dataframe
  processed_df <- process_dataframe(df)
  
  # Perform CLR transformation
  clr_tr_df_v1 <- clr_transformation(processed_df, metadata)
  clr_tr_df <- merge_by_rownames(t(clr_tr_df_v1), metadata)
  # Perform log2 transformation
  log2_tr_df <- log_transform(processed_df)
  log2_tr_df <- merge_by_rownames(log2_tr_df, metadata)
  # Perform rel. abundance
  rel_tr_df <-  rel_abundance(t(processed_df))
  rel_tr_df <- merge_by_rownames(rel_tr_df, metadata)
  # Save the result in fun_clr_list
  fun_clr_list[[df_name]] <- clr_tr_df
  fun_log2_list[[df_name]] <- log2_tr_df 
  fun_rel_list[[df_name]] <- rel_tr_df
  #plot_pcoa_and_permanova(log2_tr_df, paste0("LOG2_",df_name))
  #plot_pcoa_and_permanova(rel_tr_df, paste0("RelAbund_",df_name))
  #create_clr_pcoa_plot(t(clr_tr_df_v1), metadata, Colors, paste0("output_plots/PCOA_Plot_CLR_",df_name,".pdf"))
}
##########
##########
##########

# Loop through each dataframe in fun_rel_list
for (name1 in names(fun_rel_list)) {
  # Extract the dataframe
  #famdata <- fun_rel_list[[name]]
  
  funtable_rel <- fun_rel_list[[name1]]
  funtable_clr <- fun_clr_list[[name1]]
  funtable_log2 <- fun_log2_list[[name1]]
  
  # Create directories for each dataframe
  dir.create(name1)
  dir.create(file.path(name1, "output_tables"))
  dir.create(file.path(name1, "output_plots"))
  
  write.table(funtable_rel, file = paste(name1, "/RelAbundance_", name1, ".txt", sep = ""), sep = "\t")
  write.table(funtable_clr, file = paste(name1, "/CLRAbundance_", name1, ".txt", sep = ""), sep = "\t")
  write.table(funtable_log2, file = paste(name1, "/Log2Abundance_", name1, ".txt", sep = ""), sep = "\t")
  
  # Split status column
  status_split_funtable <- split_status_column(fun_rel_list[[name1]])
  
  plot_pcoa_and_permanova(funtable_log2, paste0(name1, "/output_plots/LOG2_",name1))
  plot_pcoa_and_permanova(funtable_rel, paste0(name1, "/output_plots/RelAbund_",name1))
  create_clr_pcoa_plot(funtable_clr[,-1], metadata, Colors, paste0(name1,"/output_plots/PCOA_Plot_CLR_",name1,".pdf"))
  
  # Put the dataframes into a list
  dataframe_list <- list(
    funtion_four_groups = fun_rel_list[[name1]],
    function_mouse_diet = status_split_funtable$df_last,
    function_human_diet = status_split_funtable$df_second_last
  )
  # Loop through the list and assign each dataframe to a new name 'famdata'
  for (name in names(dataframe_list)) {
    
    famdata <- dataframe_list[[name]]
    
    # Define the groups and their respective names
    group_vectors <- list()
    group_names <- sort(levels(as.factor(famdata$Status)))
    
    iva <- indval(famdata[,2:ncol(famdata)], famdata$Status)
    gr <- iva$maxcls[iva$pval <= 0.05]
    iv <- iva$indcls[iva$pval <= 0.05]
    pv <- iva$pval[iva$pval <= 0.05]
    fr <- apply(famdata[,2:ncol(famdata)] > 0, 2, sum)[iva$pval <= 0.05]
    indvalsummary <- data.frame(group = gr, indval = iv, pvalue = pv, freq = fr)
    indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
    subset_indvalsummary <- indvalsummary %>% filter(indval >= 0.3 & pvalue <= 0.05)
    write.table(indvalsummary, file = paste(name1, "/output_tables/labdsv_summary_", name, ".txt", sep = ""), sep = "\t")
    
    # Boruta analysis
    boruta.train <- Boruta(famdata[, 2:ncol(famdata)], as.factor(famdata$Status), pValue = 0.05, mcAdj = TRUE, maxRuns = 100, doTrace = 0, holdHistory = TRUE, getImp = getImpRfZ)
    boruta_decision <- as.data.frame(boruta.train$finalDecision)
    write.table(boruta_decision, file = paste(name1, "/output_tables/Boruta_summary_", name, ".txt", sep = ""), sep = "\t")
    
    colnames(boruta_decision)[1] <- "Decision"
    
    confirmed_df <- boruta_decision %>% filter(Decision == "Confirmed")
    famdata_t <- as.data.frame(t(famdata))
    Boruta_conf_genera <- famdata_t %>% filter(row.names(famdata_t) %in% row.names(confirmed_df))
    Boruta_conf_genera_t <- as.data.frame(t(Boruta_conf_genera))
    
    tentative_df <- boruta_decision %>% filter(Decision == "Tentative")
    Boruta_tent_genera <- famdata_t %>% filter(row.names(famdata_t) %in% row.names(tentative_df))
    Boruta_tent_genera_t <- as.data.frame(t(Boruta_tent_genera))
    
    common_rows <- intersect(rownames(subset_indvalsummary), rownames(Boruta_conf_genera))
    famdata_common_clr <- funtable_clr[, c(common_rows, "Status")]
    famdata_common_log2 <- funtable_log2[, c(common_rows, "Status")]
    famdata_common <- famdata[, c(common_rows, "Status")]
    famdata_common_clr <- famdata_common_clr %>% select(Status, everything())
    famdata_common_log2 <- famdata_common_log2 %>% select(Status, everything())
    famdata_common <- famdata_common %>% select(Status, everything())
    
    #----------ANCOM-BC
    #famdata_common will be updated by including results from ANCOM-BC
    #rm(list = c("famdata"))
    #famdata <- sp_pcoa_count
    splitted_list <- split_dataframe(famdata)
    
    se <- SummarizedExperiment(assays = list(counts = t(splitted_list$sub_count)),
                               colData = splitted_list$sub_metadata,
                               rowData = splitted_list$sub_taxa)
    
    tse <- as(se, "TreeSummarizedExperiment")
    phy <- makePhyloseqFromTreeSE(tse)
    
    out = ancombc(
      phyloseq = phy, 
      formula = "Status", 
      p_adj_method = "fdr", 
      lib_cut = 0, 
      group = "Status", 
      struc_zero = TRUE, 
      neg_lb = TRUE, 
      tol = 1e-5, 
      max_iter = 100, 
      conserve = TRUE, 
      alpha = 0.01, 
      global = TRUE
    )
    
    res <- out$res
    write.table(res, file = paste(name1, "/output_tables/ANCOM-BC_summary_", name, ".txt", sep = ""), sep = "\t")
    
    ancombc_diff_tax_list <- list()
    # Assign values to each dynamic variable based on the levels
    for (level in group_names) {
      # Create a variable name dynamically
      variable_name <- paste(level, "ancombc_diff_tax", sep = "_")
      
      # Assign values to the dynamically created variable
      if (level == group_names[1]) {
        assign(variable_name, out$res$diff_abn$taxon[out$res$diff_abn$`(Intercept)` == "TRUE"])
      } else {
        column_name <- paste("Status", level, sep = "")
        assign(variable_name, out$res$diff_abn$taxon[out$res$diff_abn[[column_name]] == "TRUE"])
      }
      ancombc_diff_tax_list[[variable_name]] <- get(variable_name)
    }
    
    # Sort the combined vector
    ancombc_diff_tax_all <- sort(unique(unlist(ancombc_diff_tax_list)))
    common_diff_taxa <- intersect(colnames(famdata_common), ancombc_diff_tax_all)
    famdata_common_v2 <- subset_dataframe(famdata_common, common_diff_taxa, c(""))
    famdata_common_clr_v2 <- subset_dataframe(famdata_common_clr, common_diff_taxa, c(""))
    famdata_common_log2_v2 <- subset_dataframe(famdata_common_log2, common_diff_taxa, c(""))
    
    famdata_common_v2 <- merge_by_rownames(famdata_common_v2, metadata)
    famdata_common_clr_v2 <- merge_by_rownames(famdata_common_clr_v2, metadata)
    famdata_common_log2_v2 <- merge_by_rownames(famdata_common_log2_v2, metadata)
    
    write.table(famdata_common_v2, file = paste0(name1, "/output_tables/diff_abund_", name1, name, ".txt"), sep = "\t")
    write.table(famdata_common_clr_v2, file = paste0(name1, "/output_tables/diff_abund_clr_", name1, name, ".txt"), sep = "\t")
    write.table(famdata_common_log2_v2, file = paste0(name1, "/output_tables/diff_abund_log2_", name1, name, ".txt"), sep = "\t")
    #edit famdata_common
    #----------------
    
    group_labdsv <- character()
    for (i in 1:length(group_names)) {
      level <- group_names[i]
      # Assign variable names dynamically
      group_vectors[[level]] <- rownames(indvalsummary)[indvalsummary$group == i]
      group_labdsv <- c(group_names)
    }
    
    # Remove all empty elements from the list
    group_vectors <- Filter(is_not_empty, group_vectors)
    subset_famdata <- function(group_vectors, famdata) {
      subsetted_data_list <- lapply(group_vectors, function(group_vector) {
        common_cols <- intersect(colnames(famdata), group_vector)
        subsetted_data <- famdata %>%
          select(Status, all_of(common_cols))
        return(subsetted_data)
      })
      return(subsetted_data_list)
    }
    # Create a PDF device to save the plots
    pdf(paste0(name1, "/output_plots/diff_group_plots_", name1, name, ".pdf"))
    
    # Loop through each vector, subset the dataframe, and generate the plot
    for (i in seq_along(group_vectors)) {
      group_vector <- group_vectors[[i]]
      group_name <- group_names[i]
      
      # Subset the dataframe using the function that handles lists
      subsetted_data_list <- subset_famdata(group_vectors, famdata_common_clr_v2)
      # Function to remove data frames with only one column
      remove_single_column_dfs <- function(list_of_dfs) {
        # Filter out data frames with more than one column
        list_of_dfs_filtered <- Filter(function(df) ncol(df) > 1, list_of_dfs)
        return(list_of_dfs_filtered)
      }
      
      # Apply the function to remove single-column data frames
      subsetted_data_list <- remove_single_column_dfs(subsetted_data_list)
      
      subsetted_data <- subsetted_data_list[[i]]
      
      create_boxplots(subsetted_data, paste0(name1, "/output_plots/", name, "_", group_name))
      
      # Reshape data for plotting
      famdata_common_long <- subsetted_data %>%
        tidyr::pivot_longer(cols = -Status, names_to = "Species", values_to = "Abundance")
      
      # Generate the plot
      p <- ggplot(famdata_common_long, aes(x = Species, y = Abundance, fill = Status)) +
        geom_boxplot(width = 0.7, position = position_dodge(width = 0.8), color = "black") +  # Adding black outline to boxplot
        geom_jitter(aes(color = Status), size = 0.5, alpha = 0.7, position = position_dodge(width = 0.8)) +
        labs(x = "Species", y = "Abundance", title = paste("Clustered Boxplot with Jittered Points -", group_name)) +
        scale_fill_manual(values = Colors) +  # Use custom colors for fill
        scale_color_manual(values = Colors) +  # Use custom colors for jittered points' outline
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        coord_flip()  # Flip axes for better readability
      
      # Print the plot to the PDF
      print(p)
    }
    dev.off()
  }
}
