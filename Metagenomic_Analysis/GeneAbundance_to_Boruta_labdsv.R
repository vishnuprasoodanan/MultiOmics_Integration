library(Boruta)
library(labdsv)
#---------------------------------------------------------------------------------------
df <- read.table("kraken2_species_abundance.txt", sep = "\t", header = T, row.names = 1)
metadata <- read.table("Metadata.txt", sep = "\t", header = T, row.names = 1)

for (col in names(df)) {
  df[[col]] <- as.numeric(df[[col]])
}

df[df <= 0.002] <- 0

zero_sum_cols <- sapply(df, function(x) sum(x) == 0)
df <- df[, !zero_sum_cols]
# Remove columns with sum zero
df <- df[, !zero_sum_cols]

# Identify rows with sum zero
zero_sum_rows <- rowSums(df) == 0

# Remove rows with sum zero
df_ref <- df[!zero_sum_rows, ]

df_ref_t <- as.data.frame(t(df_ref))

# Merge the dataframes by row names
merged_df <- merge(df_ref_t, metadata, by = "row.names", all = TRUE)

# Rename the first column to 'Name'
colnames(merged_df)[1] <- "Name"

# Set the "Name" column as rownames
rownames(merged_df) <- merged_df$Name

# Remove the "Name" column from the dataframe
merged_df <- merged_df[, !names(merged_df) %in% "Name"]

#---------Boruta Analysis
last <- as.numeric(ncol(merged_df)-1)
boruta.train <- Boruta(merged_df[,1:last], as.factor(merged_df$Status), pValue = 0.05, mcAdj = TRUE, maxRuns = 100,doTrace = 0, holdHistory = TRUE, getImp = getImpRfZ)

boruta_decision <- as.data.frame(boruta.train$finalDecision)
colnames(boruta_decision)[1] <- "Decision"

confirmed_df <- boruta_decision %>% filter(`Decision` == "Confirmed")
merged_df_t <- as.data.frame(t(merged_df))
Boruta_conf_genera <- merged_df_t %>% filter(row.names(merged_df_t) %in% row.names(confirmed_df))
Boruta_conf_genera_t <- as.data.frame(t(Boruta_conf_genera))

tentative_df <- boruta_decision %>% filter(`Decision` == "Tentative")
Boruta_tent_genera <- merged_df_t %>% filter(row.names(merged_df_t) %in% row.names(tentative_df))
Boruta_tent_genera_t <- as.data.frame(t(Boruta_tent_genera))

# Merge the dataframes by row names
Merge_Boruta_conf_genera <- merge(Boruta_conf_genera_t, metadata, by = "row.names", all = TRUE)

# Rename the first column to 'Name'
colnames(Merge_Boruta_conf_genera)[1] <- "Name"

# Set the "Name" column as rownames
rownames(Merge_Boruta_conf_genera) <- Merge_Boruta_conf_genera$Name
Merge_Boruta_conf_genera <- Merge_Boruta_conf_genera[, c("Status", setdiff(names(Merge_Boruta_conf_genera), "Status"))]
Merge_Boruta_conf_genera[, -c(1, 2)] <- sapply(Merge_Boruta_conf_genera[, -c(1, 2)], as.numeric)
# Merge the dataframes by row names
Merge_Boruta_tent_genera <- merge(Boruta_tent_genera_t, metadata, by = "row.names", all = TRUE)

# Rename the first column to 'Name'
colnames(Merge_Boruta_tent_genera)[1] <- "Name"

# Set the "Name" column as rownames
rownames(Merge_Boruta_tent_genera) <- Merge_Boruta_tent_genera$Name
Merge_Boruta_tent_genera <- Merge_Boruta_tent_genera[, c("Status", setdiff(names(Merge_Boruta_tent_genera), "Status"))]
Merge_Boruta_tent_genera[, -c(1, 2)] <- sapply(Merge_Boruta_tent_genera[, -c(1, 2)], as.numeric)

write.table(boruta.train$finalDecision, file = "P_Boruta_finaldecision.txt", sep = '\t', quote = FALSE)
write.table(boruta.train$ImpHistory, file = "P_Path_Boruta_Impotance.txt", sep = '\t', quote = FALSE)
#-------------------------------------------------------------------------------------------------------
#labdsv analysis
famdata <- merged_df
iva <- indval(famdata[,1:(ncol(famdata)-1)], famdata$Status)
gr <- iva$maxcls[iva$pval<=0.05]
iv <- iva$indcls[iva$pval<=0.05]
pv <- iva$pval[iva$pval<=0.05]
fr <- apply(famdata[,1:(ncol(famdata)-1)] >0, 2, sum)[iva$pval<=0.05]
indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
write.table (indvalsummary, file="labdsv_results.txt" , sep = "\t")
selected_genera_labdsv <- indvalsummary %>% filter(`indval` >= 0.50)
labdsv_genera <- merged_df_t %>% filter(row.names(merged_df_t) %in% row.names(selected_genera_labdsv))
# Merge the dataframes by row names
labdsv_genera_t <- as.data.frame(t(labdsv_genera))
Merge_labdsv_genera <- merge(labdsv_genera_t, metadata, by = "row.names", all = TRUE)
# Rename the first column to 'Name'
colnames(Merge_labdsv_genera)[1] <- "Name"
# Set the "Name" column as rownames
rownames(Merge_labdsv_genera) <- Merge_labdsv_genera$Name
Merge_labdsv_genera <- Merge_labdsv_genera[, c("Status", setdiff(names(Merge_labdsv_genera), "Status"))]
Merge_labdsv_genera[, -c(1, 2)] <- sapply(Merge_labdsv_genera[, -c(1, 2)], as.numeric)
#-------------------------------------------------------------------------------------------------------
#----BoxPlots
dataframes <- list(Merge_Boruta_conf_genera, Merge_Boruta_tent_genera, Merge_labdsv_genera)
# Loop through each dataframe
for (df in dataframes) {
  data <- df
  filenames <- c("Boruta_conf_genera", "Merge_Boruta_tent_genera", "Merge_labdsv_genera")
  for (name in filenames) {

    file_name1 <- paste(name, "_boxplot.pdf", sep = "_")
    pdf(file = file_name1)
    colnames(data) -> Species_name

    A_BL_Chow <- data[data$Status == "A_BL_Chow",]
    B_BL_WSD <- data[data$Status == "B_BL_WSD",]
    C_HF_Chow <- data[data$Status == "C_HF_Chow",]
    C_HF_WSD <- data[data$Status == "C_HF_WSD",]

    Colors <- c("darkolivegreen4", "red", "darkgreen", "salmon4", "wheat", "goldenrod1", "mediumpurple3", "snow4", "deepskyblue", "paleturquoise" ,"aquamarine")
    Colors1 <- c("darkolivegreen","red4", "green", "salmon", "wheat4", "goldenrod4", "mediumpurple4", "grey27", "blue", "paleturquoise4", "aquamarine4")

    plot_lst <- vector("list")
    for (i in 3:ncol(data)) {    
  
      species = Species_name[i]
      data1 = data[,c(1:2,i)]
      colnames(data1) <- c("Status", "ID", "name")
  
      A_BL_Chow_all <- A_BL_Chow[,i]
      B_BL_WSD_all <- B_BL_WSD[,i]
      C_HF_Chow_all <- C_HF_Chow[,i]
      C_HF_WSD_all <- C_HF_WSD[,i]
  
      print(species)
      species <- gsub("_", " ", species)
      P<- ggplot(data1, aes(Status, name, fill=Status))+
        ggtitle(species)+
        labs(y = "Relative-Abundance")+
        geom_boxplot(outlier.shape=NA)+ 
        geom_point(aes(colour = factor(data$Status)), position=position_jitterdodge(jitter.width = 0.5))+
        scale_color_manual(values = Colors1)+
        scale_fill_manual(values = Colors)+ theme_bw() + theme(legend.position="none")
      plot_lst[[i]] <- P
    }
    ml <- marrangeGrob(plot_lst, nrow = 1, ncol = 1)
    print(ml)
    dev.off()
  }
}
kruskalmc(data$Blautia ~ data$Status, probs=0.05)
