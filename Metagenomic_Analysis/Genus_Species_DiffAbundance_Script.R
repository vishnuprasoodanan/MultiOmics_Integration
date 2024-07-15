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

# Define parameters
Colors <- c("darkolivegreen4", "red", "darkgreen", "salmon4")
alpha_value <- 0.5
circle_alpha <- 0.8

# Function to normalize each column by dividing by the column sum
normalize_column <- function(column) {
  column_sum <- sum(column)
  return(column / column_sum)
}
#
log_transform <- function(data) {
  logcounts <- log2(data + 1)
  logcounts_matrix <- as.matrix(t(logcounts))
  return(logcounts_matrix)
}
# Function to check if a row meets the first criterion (80% of columns > 0)
check_80_percent <- function(row) {
  sum(row > 0) >= 0.8 * length(row)
}

# Function to check if a row meets the second criterion (at least one column >= 0.01)
check_at_least_one <- function(row) {
  any(row >= 0.01)
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
  
  pdf_file <- file.path("output_plots", paste0(filename, ".pdf"))
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
  
  pdf_file <- paste0("output_plots/PCOA_Plot_", output_filename, ".pdf")
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

# Function to do CLR transformatio and plot PCoA using Eucledian distance and perform PERMANOVA
create_clr_pcoa_plot <- function(counts, samples, Colors, output_file) {
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
  clr_assay <- t(assays(tse)$clr)
  clr_assay_df <- as.data.frame(clr_assay)
  
  # Calculate Euclidean distances
  euclidean_dist <- vegdist(clr_assay, method = "euclidean")
  
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
  return(clr_assay_df)
}

#------------------------------------Main code starts here----------------------------------------

# Create output directories if they don't exist
if (!dir.exists("output_tables")) dir.create("output_tables")
if (!dir.exists("output_plots")) dir.create("output_plots")


# Read data
genus_count <- read.csv("kraken2_genus_readcount.txt", sep = "\t", header = TRUE, row.names = 1)
species_count <- read.csv("kraken2_species_readcount.txt", sep = "\t", header = TRUE, row.names = 1)
metadata <- read.csv("Metadata.txt", sep = "\t", header = TRUE, row.names = 1)

# Normalize the genus and species counts
genus_abund <- as.data.frame(apply(genus_count, 2, normalize_column))
species_abund <- as.data.frame(apply(species_count, 2, normalize_column))
clr_spec_abund <- create_clr_pcoa_plot(species_count, metadata, Colors, "output_plots/CLR-transform_species.pdf")
log2_species_count <- as.data.frame(log_transform(species_count))
# Apply criteria to filter rows
rows_meeting_80_percent <- apply(genus_abund, 1, check_80_percent)
rows_meeting_at_least_one <- apply(genus_abund, 1, check_at_least_one)
filtered_genus_abund <- genus_abund[rows_meeting_80_percent & rows_meeting_at_least_one, ]
filtered_genus_abund <- filtered_genus_abund[rownames(filtered_genus_abund) != "Homo", ]

# Identify core genera
core_genera <- rownames(filtered_genus_abund)

# Create a list to store the new species abundance dataframes
spec_abund_list <- list()

# Loop through each core genus and subset species_abund
for (genus in core_genera) {
  pattern <- paste0("^[^\\w]*", genus, "\\b")
  matching_rows <- grepl(pattern, rownames(species_abund), ignore.case = FALSE)
  if (any(matching_rows)) {
    spec_abund_list[[genus]] <- species_abund[matching_rows, ]
  }
}

# Create new dataframes by transposing and matching with metadata
new_dataframes_list <- list()

for (genus in names(spec_abund_list)) {
  transposed_df <- t(spec_abund_list[[genus]])
  transposed_df <- as.data.frame(transposed_df)
  transposed_df$sample_id <- rownames(transposed_df)
  
  matching_rows <- intersect(rownames(metadata), transposed_df$sample_id)
  if (length(matching_rows) > 0) {
    merged_df <- merge(metadata[matching_rows, 1, drop = FALSE], transposed_df, by.x = "row.names", by.y = "sample_id")
    rownames(merged_df) <- merged_df$Row.names
    merged_df <- merged_df[ , -1]  # Remove the Row.names column
    new_dataframes_list[[genus]] <- merged_df
  }
}

# Modify elements in new_dataframes_list to create two_group_df_list
two_group_df_list <- lapply(new_dataframes_list, function(df) {
  if ("Status" %in% colnames(df)) {
    df$Status <- sapply(strsplit(as.character(df$Status), "_"), function(x) tail(x, 1))
  }
  return(df)
})

# Perform Wilcoxon rank sum test and identify significant rownames
significant_rownames <- c()

for (df in two_group_df_list) {
  for (gene in colnames(df)[-1]) {  # Exclude the Status column
    wilcox_test <- wilcox.test(df[[gene]] ~ df$Status)
    if (!is.na(wilcox_test$p.value) && wilcox_test$p.value <= 0.05) {
      significant_rownames <- c(significant_rownames, gene)
    }
  }
}

# Remove duplicates
significant_rownames <- unique(significant_rownames)

# Subset each dataframe in two_group_df_list by selecting significant rownames
subsetted_dataframes_list <- lapply(two_group_df_list, function(df) {
  valid_columns <- intersect(c("Status", significant_rownames), colnames(df))
  return(df[, valid_columns, drop = FALSE])
})

# Calculate average value of abundance for each group in 'Status' for each subsetted dataframe
averages <- lapply(subsetted_dataframes_list, function(df) {
  status_groups <- unique(df$Status)
  avg_values <- data.frame(Status = status_groups)
  
  high_chow <- character(0)
  high_wsd <- character(0)
  
  for (col in colnames(df)[-1]) {  # Exclude the Status column
    diff_value <- mean(df[df$Status == 'Chow', col], na.rm = TRUE) - mean(df[df$Status == 'WSD', col], na.rm = TRUE)
    if (diff_value > 0) {
      high_chow <- c(high_chow, col)
    } else if (diff_value < 0) {
      high_wsd <- c(high_wsd, col)
    }
    
    means <- sapply(status_groups, function(group) {
      mean(df[df$Status == group, col], na.rm = TRUE)
    })
    avg_values[[col]] <- means
  }
  
  return(list(avg_values = avg_values, high_chow = high_chow, high_wsd = high_wsd))
})

# Create lists to store subsets
high_chow_df_list <- list()
high_wsd_df_list <- list()

# Loop through each dataframe in subsetted_dataframes_list
for (i in seq_along(subsetted_dataframes_list)) {
  df <- subsetted_dataframes_list[[i]]
  
  high_chow <- averages[[i]]$high_chow
  high_wsd <- averages[[i]]$high_wsd
  
  high_chow_df <- new_dataframes_list[[i]][, c("Status", high_chow), drop = FALSE]
  high_wsd_df <- new_dataframes_list[[i]][, c("Status", high_wsd), drop = FALSE]
  
  high_chow_df_list[[i]] <- high_chow_df
  high_wsd_df_list[[i]] <- high_wsd_df
  
  names(high_chow_df_list)[i] <- names(subsetted_dataframes_list)[i]
  names(high_wsd_df_list)[i] <- names(subsetted_dataframes_list)[i]
}

# List of all dataframes
all_df_list <- list(
  high_chow = high_chow_df_list,
  high_wsd = high_wsd_df_list
)

# Loop through each dataframe in all_df_list and create boxplots
for (list_name in names(all_df_list)) {
  for (i in seq_along(all_df_list[[list_name]])) {
    df <- all_df_list[[list_name]][[i]]
    df_name <- names(all_df_list[[list_name]])[i]
    
    create_boxplots(df, paste0(list_name, "_", df_name))
  }
}

# Save filtered_genus_abund to the output_tables directory
write.table(filtered_genus_abund, file = "output_tables/core_genus_abund.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# Transpose the species_abund dataframe and replace values < 0.002 with zero
transposed_species_abund <- t(species_abund)
transposed_species_count <- t(species_count)
transposed_species_count[transposed_species_count < 100] <- 0
transposed_species_count <- transposed_species_count[, colSums(transposed_species_count == 0) <= 4]
columns_sel <- colnames(transposed_species_count)

sp_pcoa_abund <- as.data.frame(subset_dataframe(transposed_species_abund, columns_sel, c("")))
sp_pcoa_abund <- merge_by_rownames(sp_pcoa_abund, metadata)
# Print dimensions of the final dataframe for verification
print(dim(sp_pcoa_abund))
print(head(sp_pcoa_abund))
plot_pcoa_and_permanova(sp_pcoa_abund, "sp_pcoa_abund")
# Perform PCoA and PERMANOVA for each dataframe in new_dataframes_list
for (i in seq_along(new_dataframes_list)) {
  df <- new_dataframes_list[[i]]
  df_name <- names(new_dataframes_list)[i]
  
  if (ncol(df) >= 5) {
    plot_pcoa_and_permanova(df, df_name)
  } else {
    cat("Skipping dataframe", df_name, "because it has less than 5 columns.\n")
  }
}


sp_pcoa_count <- as.data.frame(t(subset_dataframe(species_count ,c(""), colnames(sp_pcoa_abund))))
sp_pcoa_count <- merge_by_rownames(sp_pcoa_count, metadata)
write.table(sp_pcoa_count, file = "output_tables/filtered_species_count.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(sp_pcoa_abund, file = "output_tables/filtered_species_abund.txt", sep = "\t", quote = FALSE, row.names = TRUE)

#----------------------------CLR transformation
selected_rows <- intersect(colnames(sp_pcoa_abund), rownames(species_count))
species_count_subset <- subset_dataframe(species_count, c(""), selected_rows)
clr_core_spec <- create_clr_pcoa_plot(species_count_subset, metadata, Colors, "output_plots/CLR-transform_subset_species.pdf")

# Subset clr_spec_abund by matching column names
selected_columns <- intersect(colnames(sp_pcoa_abund), colnames(clr_spec_abund))
clr_spec_abund_subset <- subset_dataframe(clr_spec_abund, selected_columns, c(""))
clr_spec_abund_subset <- merge_by_rownames(clr_spec_abund_subset, metadata)

write.table(clr_spec_abund_subset, file = "output_tables/filtered_species_clr.txt", sep = "\t", quote = FALSE, row.names = TRUE)

#----------------------------Log2 transformation

# Apply log2 transformation to all values in the dataframe
log2_spec_abund_subset <- subset_dataframe(log2_species_count, selected_columns, c(""))
log2_spec_abund_subset <- merge_by_rownames(log2_spec_abund_subset, metadata)
write.table(log2_spec_abund_subset, file = "output_tables/filtered_species_log2.txt", sep = "\t", quote = FALSE, row.names = TRUE)

#------------------------------differential abundance analysis
# labdsv analysis

famdata <- sp_pcoa_abund
iva <- indval(famdata[,2:ncol(famdata)], famdata$Status)
gr <- iva$maxcls[iva$pval <= 0.05]
iv <- iva$indcls[iva$pval <= 0.05]
pv <- iva$pval[iva$pval <= 0.05]
fr <- apply(famdata[,2:ncol(famdata)] > 0, 2, sum)[iva$pval <= 0.05]
indvalsummary <- data.frame(group = gr, indval = iv, pvalue = pv, freq = fr)
indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
subset_indvalsummary <- indvalsummary %>% filter(indval >= 0.5 & pvalue <= 0.05)
write.table(indvalsummary, file = "output_tables/labdsv_results.txt", sep = "\t")

# Boruta analysis
boruta.train <- Boruta(famdata[, 2:ncol(famdata)], as.factor(famdata$Status), pValue = 0.05, mcAdj = TRUE, maxRuns = 100, doTrace = 0, holdHistory = TRUE, getImp = getImpRfZ)
boruta_decision <- as.data.frame(boruta.train$finalDecision)
write.table(boruta_decision, file = "output_tables/boruta_results.txt", sep = "\t")

colnames(boruta_decision)[1] <- "Decision"

confirmed_df <- boruta_decision %>% filter(Decision == "Confirmed")
famdata_t <- as.data.frame(t(famdata))
Boruta_conf_genera <- famdata_t %>% filter(row.names(famdata_t) %in% row.names(confirmed_df))
Boruta_conf_genera_t <- as.data.frame(t(Boruta_conf_genera))

tentative_df <- boruta_decision %>% filter(Decision == "Tentative")
Boruta_tent_genera <- famdata_t %>% filter(row.names(famdata_t) %in% row.names(tentative_df))
Boruta_tent_genera_t <- as.data.frame(t(Boruta_tent_genera))

common_rows <- intersect(rownames(subset_indvalsummary), rownames(Boruta_conf_genera))
famdata_common <- clr_spec_abund_subset[, c(common_rows, "Status")]
famdata_common <- famdata_common %>% select(Status, everything())

#----------ANCOM-BC
#famdata_common will be updated by including results from ANCOM-BC
rm(list = c("famdata"))
famdata <- sp_pcoa_count
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
write.table(res, file = "output_tables/ancom-bc_results.txt", sep = "\t")
# Initialize vectors
A_BL_Chow_ancombc_diff_tax <- c()
B_BL_WSD_ancombc_diff_tax <- c()
C_HF_Chow_ancombc_diff_tax <- c()
D_HF_WSD_ancombc_diff_tax <- c()

# Assign values to A_BL_Chow_ancombc_diff_tax if the second column is "TRUE"
A_BL_Chow_ancombc_diff_tax <- out$res$diff_abn$taxon[out$res$diff_abn$`(Intercept)` == "TRUE"]
B_BL_WSD_ancombc_diff_tax <- out$res$diff_abn$taxon[out$res$diff_abn$StatusB_BL_WSD == "TRUE"]
C_HF_Chow_ancombc_diff_tax <- out$res$diff_abn$taxon[out$res$diff_abn$StatusC_HF_Chow == "TRUE"]
D_HF_WSD_ancombc_diff_tax <- out$res$diff_abn$taxon[out$res$diff_abn$StatusD_HF_WSD == "TRUE"]

# Combine all vectors into a single vector and ensure uniqueness
ancombc_diff_tax_all <- unique(c(A_BL_Chow_ancombc_diff_tax, 
                                 B_BL_WSD_ancombc_diff_tax, 
                                 C_HF_Chow_ancombc_diff_tax, 
                                 D_HF_WSD_ancombc_diff_tax))

# Sort the combined vector
ancombc_diff_tax_all <- sort(ancombc_diff_tax_all)
common_diff_taxa <- intersect(colnames(famdata_common), ancombc_diff_tax_all)
famdata_common_v2 <- subset_dataframe(famdata_common, common_diff_taxa, c(""))
famdata_common_v2 <- merge_by_rownames(famdata_common_v2, metadata)
#edit famdata_common
#----------------
A_BL_Chow_group <- rownames(indvalsummary)[indvalsummary$group == 1]
B_BL_WSD_group <- rownames(indvalsummary)[indvalsummary$group == 2]
C_HF_Chow_group <- rownames(indvalsummary)[indvalsummary$group == 3]
D_HF_WSD_group <- rownames(indvalsummary)[indvalsummary$group == 4]

# Function to subset famdata_common by common columns with group vector
subset_famdata <- function(group_vector, famdata) {
  common_cols <- intersect(colnames(famdata), group_vector)
  subsetted_data <- famdata %>%
    select(Status, all_of(common_cols))
  return(subsetted_data)
}

# Define the groups and their respective names
group_vectors <- list(A_BL_Chow_group, B_BL_WSD_group, C_HF_Chow_group, D_HF_WSD_group)
group_names <- c("A_BL_Chow_group", "B_BL_WSD_group", "C_HF_Chow_group", "D_HF_WSD_group")

# Create a PDF device to save the plots
pdf("output_plots/group_plots.pdf")

# Loop through each vector, subset the dataframe, and generate the plot
for (i in seq_along(group_vectors)) {
  group_vector <- group_vectors[[i]]
  group_name <- group_names[i]
  
  # Subset the dataframe
  subsetted_data <- subset_famdata(group_vector, famdata_common_v2)
  create_boxplots(subsetted_data, paste0(group_name, "_", group_name))
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
# Close the PDF device
dev.off()

#---------------------------Perform Mantel test

log2_spec_abund_subset 
clr_spec_abund_subset
sp_pcoa_abund

clr_sp_euclidean_dist <- vegdist(clr_spec_abund_subset[,-1], method = "euclidean")
log2_bray_distances <- vegdist(log2_spec_abund_subset[, -1], method = "bray")
rel_bray_distances <- vegdist(sp_pcoa_abund[, -1], method = "bray")

# List of dataframes (excluding the first column of each)
ab_df_list <- list(
  clr_spec_abund_subset = clr_spec_abund_subset[,-1],
  log2_spec_abund_subset = log2_spec_abund_subset[,-1],
  rel_spec_abund_subset = sp_pcoa_abund[,-1]
)

dist_df_list <- list(
  clr_sp_euclidean_dist = clr_sp_euclidean_dist,
  log2_bray_distances = log2_bray_distances,
  rel_bray_distances = rel_bray_distances
)

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

# Perform Mantel tests and store results
results_df <- perform_mantel_test(ab_df_list)
results_dist <- perform_mantel_test(dist_df_list)

print(results_df)
print(results_dist)

write.table(results_df, file = "output_tables/abundance_mantel_test.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(results_dist, file = "output_tables/distance_mantel_test.txt", sep = "\t", quote = FALSE, row.names = TRUE)
