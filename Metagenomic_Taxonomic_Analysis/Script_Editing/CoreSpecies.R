KO <- read.table(file = "species_abundance_nnumeric.txt", sep = '\t', colClasses = "numeric") # data without row and column names
KO[KO < 0.002] <- 0
KO <- KO[rowSums(KO == 0) <= 6,]

KO_test <- read.table(file = "species_abundance_all_data.txt", sep = '\t', header = TRUE, row.names = 1)
KO_test[KO_test < 0.002] <- 0
KO_test <- KO_test[rowSums(KO_test == 0) <= 6,]
write.table(KO_test, file = "Core_Species_80.txt", quote = FALSE, sep = '\t', row.names = TRUE)

Blautia_rows <- grepl("Blautia", row.names(KO_test))
# subset the dataframe to include only rows that contain "abcxy"
Blautia_df <- KO_test[Blautia_rows,]
# print the resulting dataframe
write.table(Blautia_df, file = "Blautia_species_abundance.txt", quote = FALSE, sep = '\t', row.names = TRUE)

Bacteroides_rows <- grepl("Bacteroides", row.names(KO_test))
# subset the dataframe to include only rows that contain "abcxy"
Bacteroides_df <- KO_test[Bacteroides_rows,]
# print the resulting dataframe
write.table(Bacteroides_df, file = "Bacteroides_species_abundance.txt", quote = FALSE, sep = '\t', row.names = TRUE)

Erysipelatoclostridium_rows <- grepl("Erysipelatoclostridium", row.names(KO_test))
# subset the dataframe to include only rows that contain "abcxy"
Erysipelatoclostridium_df <- KO_test[Erysipelatoclostridium_rows,]
# print the resulting dataframe
write.table(Erysipelatoclostridium_df, file = "Erysipelatoclostridium_species_abundance.txt", quote = FALSE, sep = '\t', row.names = TRUE)

Parabacteroides_rows <- grepl("Parabacteroides", row.names(KO_test))
# subset the dataframe to include only rows that contain "abcxy"
Parabacteroides_df <- KO_test[Parabacteroides_rows,]
# print the resulting dataframe
write.table(Parabacteroides_df, file = "Parabacteroides_species_abundance.txt", quote = FALSE, sep = '\t', row.names = TRUE)

Lachnocostridium_rows <- grepl("Lachnocostridium", row.names(KO_test))
# subset the dataframe to include only rows that contain "abcxy"
Lachnocostridium_df <- KO_test[Lachnocostridium_rows,]
# print the resulting dataframe
write.table(Lachnocostridium_df, file = "Lachnocostridium_species_abundance.txt", quote = FALSE, sep = '\t', row.names = TRUE)

KO_test_T <- as.data.frame(t(KO_test))
KO_test_1 <- as.data.frame(KO_test_T[2:367])
rownames(KO) <- colnames(KO_test_1)
colnames(KO) <- row.names(KO_test_1)

KO_proportions <- KO
KO_proportions1 <- as.data.frame(t(KO_proportions))

class <- KO_test_T$Status
Bray_pcoa <-pcoa(vegdist(KO_proportions1, "bray"))
Bray_pcoa$values[1:2,]
mds.var.per = round(Bray_pcoa$values$Eigenvalues/sum(Bray_pcoa$values$Eigenvalues)*100, 1)
Bray_PCoA_MATRIX <- Bray_pcoa$vectors[,1:2]
Bray_PCoA_MATRIX <- as.data.frame(Bray_PCoA_MATRIX)

Bray_distances <-vegdist(KO_proportions1, "bray")
adonis2(Bray_distances ~ KO_test_T$Status)

Bray_PCoA_MATRIX_New <- cbind(Bray_PCoA_MATRIX, class)
write.table(Bray_PCoA_MATRIX_New, file = "PCA_data", quote = FALSE, sep = '\t')

pc <- c(1,2)
jpeg("PCOA_count_data_v2.jpg", height = 10, width = 10, units = 'in', res = 600)

plot(Bray_pcoa$vectors[,1:2], bg=c("darkolivegreen4", "red", "darkgreen", "salmon4")[as.factor(Bray_PCoA_MATRIX_New$class)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
ordiellipse(Bray_pcoa$vectors[,1:2], Bray_PCoA_MATRIX_New$class, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 70, col = c("darkolivegreen4", "red", "darkgreen", "salmon4"))
#ordispider(Bray_pcoa$vectors[,1:2], Bray_PCoA_MATRIX_New$class, lty=3, spider ="centroid", lwd=1, col="black")
legend("topright", legend = c("A_BL_Chow", "B_BL_WSD", "C_HF_Chow", "C_HF_WSD"), col = c("darkolivegreen4", "red", "darkgreen", "salmon4"),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
text(Bray_pcoa$vectors[,1:2], labels=as.factor(rownames(KO_proportions1)), cex=0.6, font=1, pos=1)
#abline(h=0, v=0, col = "gray60")
dev.off ()
