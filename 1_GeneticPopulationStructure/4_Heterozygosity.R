library(poolfstat)

# Define pool names (EXACT order from sync file)
poolnames <- c("CBA", "CGB", "CPA", "HBA", "HGB", "HPA",
               "CBB", "CHA", "CPB", "HBB", "HHA", "HPB",
               "CGA", "CHB", "G1", "HGA", "HHB")

# Define corresponding pool sizes (EXACT order from sync file)
poolsizes <- c(60, 56, 60, 54, 54, 60,
               44, 60, 60, 60, 46, 54,
               48, 60, 44, 60, 54)


# Convert sync file to PoolFstat format
pf_data <- popsync2pooldata(
  sync.file = "../../0_data/SNPs_noRepeats.sync",  # Your sync file
  poolsizes = poolsizes,              # Your exact pool sizes
  poolnames = poolnames,              # Your exact pool names
  min.rc = 1,                          # Minimum read count per SNP
  min.cov.per.pool = 15,                # Minimum coverage per pool
  max.cov.per.pool = 500,              # Maximum coverage per pool
  min.maf = 0.05,                      # Minimum minor allele frequency
  noindel = TRUE,                      # Ignore indels
  nlines.per.readblock = 1e+06,        # Read large blocks to optimize speed
  nthreads = 4                         # Use 4 threads for faster processing
)

# Print summary of the dataset
summary(pf_data)

# calculate fstats
res.fstats = compute.fstats(pf_data)



#############Heterozygosity
# visualize heterozygosity
# Load required libraries
library(ggplot2)
library(dplyr)

# Convert to data frame and sort the orders
heterozygosity_df <- data.frame(Population = poolnames, H = res.fstats@heterozygosities)
sorted_heterozygosity_df <- heterozygosity_df %>% arrange(Population)  # Ascending order

# Define color scheme
cold_colors <- colorRampPalette(c("#c6dbef", "#2171b5"))(4)  # Blue shades
hot_colors <- colorRampPalette(c("#fcae91", "#cb181d"))(4)   # Red shades
g1_color <- "grey50"  # Grey for G1

# Assign colors by repeating each treatment color for both replicates (A/B)
cold_treatments <- c("CBA", "CBB", "CPA", "CPB", "CGA", "CGB", "CHA", "CHB")
hot_treatments  <- c("HBA", "HBB", "HPA", "HPB", "HGA", "HGB", "HHA", "HHB")

# Repeat colors for replicates
cold_treatment_colors <- rep(cold_colors, each=2)  # Ensures replicates share same shade
hot_treatment_colors  <- rep(hot_colors, each=2)

# Assign colors to treatments
treatment_colors <- c(setNames(cold_treatment_colors, cold_treatments),
                      "G1" = g1_color,
                      setNames(hot_treatment_colors, hot_treatments))  # Add G1 color

# Ensure the Population column is treated as a factor with correct order
sorted_heterozygosity_df$Population <- factor(sorted_heterozygosity_df$Population, levels = names(treatment_colors))

# Plot
ggplot(sorted_heterozygosity_df, aes(x = Population, y = Estimate, fill = Population)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = treatment_colors) +
  theme_minimal() +
  labs(title = "Heterozygosity Across Treatments", 
       x = "Pools", y = "Heterozygosity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none",
        panel.grid.major.x = element_blank())  # Remove vertical grid lines



#############pairwise Fst heatmap
# Load required libraries
library(ggplot2)
library(reshape2)
library(pheatmap)

# Extract pairwise FST matrix
fst_matrix <- res.fstats@pairwise.fst


# Define the fixed order of populations (same as heterozygosity figure)
pop_order <- c("CBA", "CBB", "CPA", "CPB", "CGA", "CGB", "CHA", "CHB", 
               "G1", 
               "HBA", "HBB", "HPA", "HPB", "HGA", "HGB", "HHA", "HHB")

# Reorder matrix according to the defined order
fst_matrix <- fst_matrix[pop_order, pop_order]

write.table(fst_matrix,file="fst.txt", quote=FALSE, sep="\t") 

# Remove upper triangle to keep only the lower triangular matrix
fst_matrix[upper.tri(fst_matrix)] <- NA

# Assign colors by repeating each treatment color for both replicates (A/B)
cold_treatments <- c("CBA", "CBB", "CPA", "CPB", "CGA", "CGB", "CHA", "CHB")
hot_treatments  <- c("HBA", "HBB", "HPA", "HPB", "HGA", "HGB", "HHA", "HHB")
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

# Create annotation data frame for row and column colors
population_annotations <- data.frame(Treatment = pop_order)
rownames(population_annotations) <- pop_order


# Generate FST heatmap (NO CLUSTERING)
pheatmap(fst_matrix, 
         color = colorRampPalette(c("yellow", "orange", "red"))(50),  # FST color scale
         annotation_row = population_annotations,  
         annotation_col = population_annotations,
         annotation_colors = list(Treatment = annotation_colors),  
         cluster_rows = FALSE,  # No clustering
         cluster_cols = FALSE,  # No clustering
         display_numbers = FALSE,  
         border_color = NA,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 10)
