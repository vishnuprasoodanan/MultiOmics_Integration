library(MintTea)
library(openxlsx)

#data('test_data')
minttea_input<- read.csv("Mint_tea_input.txt", sep = "\t", header = TRUE)
metadata <- read.csv("ModifiedMetadata.txt", sep = "\t", header = TRUE, row.names = 1)

# Loop through each column of 'metadata'
for (col in colnames(metadata)) {
  # Save the column name in 'group_name'
  group_name <- col
  
  # Replace the last column of 'minttea_input' with this column
  minttea_input[, ncol(minttea_input)] <- metadata[[group_name]]
  
  # Rename the last column to 'disease_state'
  colnames(minttea_input)[ncol(minttea_input)] <- "disease_state"
  write.table (minttea_input, file = paste0(group_name,"_minttea_input.txt"), sep = "\t")
  minttea_results <- MintTea(minttea_input, view_prefixes = c('P', 'M', 'E', 'C'))
  #output_file <- "Minttea_BL_Chow.txt"
  
  # Initialize maximum values for a and b
  max_a <- 0
  max_b <- 0
  
  # Find the maximum values of 'a' and 'b'
  for (a in seq_along(minttea_results)) {
    if (is.list(minttea_results[[a]])) {
      max_a <- max(max_a, a)
      for (b in seq_along(minttea_results[[a]])) {
        max_b <- max(max_b, b)
      }
    }
  }
  
  # Print maximum values of 'a' and 'b'
  cat("Maximum value of 'a':", max_a, "\n")
  cat("Maximum value of 'b':", max_b, "\n")
  
  # Create a new workbook
  wb <- createWorkbook()
  
  # Loop through 'a' and 'b' and print outputs in separate sheets
  for (a in seq_len(max_a)) {
    if (is.list(minttea_results[[a]])) {
      for (b in seq_len(max_b)) {
        if (!is.null(minttea_results[[a]][[b]])) {
          sheet_name <- paste0(group_name,"a", a, "_b", b)
          addWorksheet(wb, sheet_name)
          
          # Capture the formatted print output
          formatted_output <- capture.output(print(minttea_results[[a]][[b]]))
          
          # Write the formatted output to the sheet
          writeData(wb, sheet_name, formatted_output, startCol = 1, startRow = 1, colNames = FALSE)
          
          cat("Printed results of minttea_results[[", a, "]][[", b, "]] to sheet", sheet_name, "\n")
        }
      }
    }
  }
  
  # Save the workbook
  saveWorkbook(wb, paste0(group_name, "_minttea_results.xlsx"), overwrite = TRUE)
  
  # Print message to indicate the task is complete
  cat("All results have been written to minttea_results.xlsx\n")
}

