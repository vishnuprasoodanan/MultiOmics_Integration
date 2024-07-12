# Load required libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(labdsv)
library(RColorBrewer)
library(vegan)
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

# Function to create boxplots for each dataframe and save as PDF
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
        
        empty_plots <- 9 - length(plots_subset)
        for (i in 1:empty_plots) {
          plots_subset <- c(plots_subset, NULL)
        }
      } else {
        num_cols <- 3
        num_rows <- 3
      }
      
      plot_grid <- do.call(grid.arrange, c(plots_subset, ncol = num_cols))
      print(plot_grid)
    }
  }
  
  dev.off()
  cat(paste("PDF file", pdf_file, "created successfully.\n"))
}

# Function to plot PCoA and perform PERMANOVA
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
  
  output_data_file <- paste0("output_tables/CAZY_Healthy1lk_SE_Spainremoved_PCA_data_", output_filename, ".txt")
  write.table(bray_pcoa_matrix_new, file = output_data_file, quote = FALSE, sep = '\t')
  
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

#----------------------------------------------------Main code starts here----------------------------------------------------

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
write.table(filtered_genus_abund, file = "output_tables/filtered_genus_abund.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# Transpose the species_abund dataframe and replace values < 0.002 with zero
transposed_species_abund <- t(species_abund)
transposed_species_abund[transposed_species_abund < 0.00001] <- 0

# Subset transposed_species_abund by excluding columns with more than 4 zeros
subsetted_species_abund <- transposed_species_abund[, colSums(transposed_species_abund == 0) <= 4]

# Match rownames of transposed_species_abund with rownames of metadata and add 'Status' column
matching_rows <- intersect(rownames(metadata), rownames(subsetted_species_abund))
sp_pcoa_abund <- as.data.frame(subsetted_species_abund[matching_rows, ])
sp_pcoa_abund$Status <- metadata[matching_rows, "Status"]
sp_pcoa_abund <- move_status_to_first(sp_pcoa_abund)

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

#------------------------------------------------differential abundance analysis-------------------------------------------------
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
write.table(indvalsummary, file = "labdsv_results.txt", sep = "\t")

# Boruta analysis
boruta.train <- Boruta(famdata[, 2:ncol(famdata)], as.factor(famdata$Status), pValue = 0.05, mcAdj = TRUE, maxRuns = 100, doTrace = 0, holdHistory = TRUE, getImp = getImpRfZ)
boruta_decision <- as.data.frame(boruta.train$finalDecision)
colnames(boruta_decision)[1] <- "Decision"

confirmed_df <- boruta_decision %>% filter(Decision == "Confirmed")
famdata_t <- as.data.frame(t(famdata))
Boruta_conf_genera <- famdata_t %>% filter(row.names(famdata_t) %in% row.names(confirmed_df))
Boruta_conf_genera_t <- as.data.frame(t(Boruta_conf_genera))

tentative_df <- boruta_decision %>% filter(Decision == "Tentative")
Boruta_tent_genera <- famdata_t %>% filter(row.names(famdata_t) %in% row.names(tentative_df))
Boruta_tent_genera_t <- as.data.frame(t(Boruta_tent_genera))

common_rows <- intersect(rownames(subset_indvalsummary), rownames(Boruta_conf_genera))
famdata_common <- famdata[, c(common_rows, "Status")]
famdata_common <- famdata_common %>% select(Status, everything())

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
  subsetted_data <- subset_famdata(group_vector, famdata_common)
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
