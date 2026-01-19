# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read allele frequency file
allele_freq <- read.table("allele_frequencies_and_coverage.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")

# Check the structure
head(allele_freq)

# Select only the allele frequency columns (drop CHROM, POS, REF, ALT)
AF_matrix <- allele_freq %>%
  select(starts_with("AF_")) %>%  # Select only allele frequency columns
  na.omit()  # Remove rows with missing values

# Transpose the matrix so pools become rows (samples) and SNPs become columns (features)
AF_matrix_t <- t(AF_matrix)

# Set pool names as row names
rownames(AF_matrix_t) <- gsub("AF_", "", colnames(AF_matrix))  # Remove "AF_" prefix
rownames(AF_matrix_t) <- gsub("pool_", "", rownames(AF_matrix_t))  # Remove "pool_" prefix

# Convert to numeric
AF_matrix_t <- as.data.frame(AF_matrix_t)
AF_matrix_t[] <- lapply(AF_matrix_t, as.numeric)  # Ensure all values are numeric

# Check structure
dim(AF_matrix_t)  # Should be (17 pools × 1787629 SNPs)

# Run PCA
pca_results <- prcomp(AF_matrix_t, center=TRUE, scale=TRUE)

# Get PCA scores
pca_df <- as.data.frame(pca_results$x)

# Add pool names as a new column
pca_df$Pool <- rownames(AF_matrix_t)

# Check PCA summary
summary(pca_results)


# Assign colors by repeating each treatment color for both replicates (A/B)
cold_treatments <- c("CBA", "CBB", "CGA", "CGB", "CHA", "CHB", "CPA", "CPB")
hot_treatments  <- c("HBA", "HBB", "HGA", "HGB", "HHA", "HHB", "HPA", "HPB")
g1_treatment    <- "G1"

# Define color scheme
cold_colors <- colorRampPalette(c("#c6dbef", "#2171b5"))(4)  # Blue shades
hot_colors <- colorRampPalette(c("#fcae91", "#cb181d"))(4)   # Red shades
g1_color <- "grey50"  # Grey for G1


# Repeat colors for replicates
cold_treatment_colors <- rep(cold_colors, each=2)  # Ensures replicates share same shade
hot_treatment_colors  <- rep(hot_colors, each=2)

# Assign colors to treatments
annotation_colors <- c(setNames(cold_treatment_colors, cold_treatments),
                       "G1" = g1_color,
                       setNames(hot_treatment_colors, hot_treatments))  # Add G1 color


# Match pools with colors
pca_df$Color <- annotation_colors[pca_df$Pool]


# Load necessary libraries
library(ggplot2)

# Create PCA plot with labels
ggplot(pca_df, aes(x=PC1, y=PC2, color=Pool, label=Pool)) +
  geom_point(size=6, alpha=0.7) +  # Plot points
  geom_text(aes(label=Pool), vjust=-1, hjust=0.5, size=4) +  # Add text annotations above points
  scale_color_manual(values=annotation_colors) +  # Use predefined colors
  labs(title="PCA of Allele Frequencies",
       x=paste0("PC1 (", round(summary(pca_results)$importance[2,1] * 100, 1), "%)"),
       y=paste0("PC2 (", round(summary(pca_results)$importance[2,2] * 100, 1), "%)")) +
  theme_minimal() +
  theme(legend.position="none")  # Hide legend since labels are in the plot
