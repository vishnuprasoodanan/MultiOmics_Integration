pdf(file = "Cluster90_Final_ScatterPlot.pdf", height = 5, width = 10)
data <- read.table(file = "Cluster90_represetatives.txt", sep = "\t", header = T, row.names = 1)
# Add a new column "Size" based on the conditions
data$status <- ifelse(data$completeness >= 95 & data$contamination <= 5, "HighlyComplete",
                      ifelse(data$completeness >= 90 & data$contamination <= 10, "NearlyComplete",
                             ifelse(data$completeness >= 70 & data$completeness < 90 & data$contamination <= 10, "MediumQuality",
                                    ifelse(data$completeness >= 50 & data$completeness < 70 & data$contamination <= 10, "PartialBins",NA))))

# Change point shapes and colors
ggplot(data, aes(x=completeness, y=contamination, color=status)) +
  geom_point(alpha=0.5, size=4, shape=16)+
  scale_color_manual(values=c("red", "darkolivegreen4", "salmon4", "deepskyblue")) +
  theme_bw()

dev.off ()
