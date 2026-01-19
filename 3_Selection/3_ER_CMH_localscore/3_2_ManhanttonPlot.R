library(data.table)
library(ggplot2)
library(dplyr)

# Define treatments and pollinator/temperature labels
treatments <- c("CB", "HB", "CP", "HP", "CG", "HG", "CH", "HH")
pollinator_labels <- c("B" = "Bee", "P" = "Butterfly", "G" = "Both", "H" = "Hand")
temperature_labels <- c("C" = "Cold", "H" = "Hot")

# Read and combine all filtered data
all_data <- rbindlist(lapply(treatments, function(tr) {
  df <- fread(paste0("cmh_", tr, ".txt"))  # Already pre-filtered to first 3 columns
  setnames(df, c("CHROM", "POS", "pval", "qval"))
  df$treatment <- tr
  df$temp <- substr(tr, 1, 1)
  df$pollinator <- substr(tr, 2, 2)
  return(df)
}))
all_data$pval[all_data$pval==0]=1e-16

# Ensure chromosome order as A01–A10 (factor)
chrom_levels <- paste0("A", sprintf("%02d", 1:10))
all_data <- all_data %>%
  mutate(
    logp = -log10(pval),
    CHROM = factor(CHROM, levels = chrom_levels)
  )

# Prepare cumulative positions (offset) for chromosomes
chrom_lengths <- all_data %>%
  group_by(CHROM) %>%
  summarise(max_pos = max(POS)) %>%
  mutate(offset = cumsum(lag(max_pos, default = 0)))

# Add global position for plotting
all_data <- all_data %>%
  left_join(chrom_lengths, by = "CHROM") %>%
  mutate(global_pos = POS + offset)

# Define chromosome colors (alternating black/grey)
chrom_colors <- rep(c("black", "grey"), length.out = length(chrom_levels))
names(chrom_colors) <- chrom_levels

###PLOT THRESHOLD FDR<0.05, remove the extreme
# Plot
# Set clipping threshold
clip_threshold <- 60

# Add a column to flag capped points
all_data <- all_data %>%
  mutate(
    logp_clipped = pmin(logp, clip_threshold),
    clipped = logp > clip_threshold
  )

# Adjust plotting
q <- ggplot(all_data, aes(x = global_pos, y = logp_clipped)) +
  # Background points (not significant)
  geom_point(
    data = subset(all_data, qval >= 0.05),
    aes(color = CHROM),
    size = 0.5
  ) +
  scale_color_manual(values = chrom_colors, guide = "none") +
  
  # Significant points (FDR < 0.05)
  geom_point(
    data = subset(all_data, qval < 0.05),
    aes(x = global_pos, y = logp_clipped),
    color = "red", size = 0.8
  ) +
  
  # Add special marker for capped points (optional, e.g., triangles for clipped)
  geom_point(
    data = subset(all_data, clipped),
    aes(x = global_pos, y = clip_threshold),
    shape = 17,  # Triangle shape
    color = "blue", size = 1.5
  ) +
  
  scale_x_continuous(
    breaks = chrom_lengths$offset + chrom_lengths$max_pos / 2,
    labels = chrom_lengths$CHROM
  ) +
  facet_grid(pollinator_labels[pollinator] ~ temperature_labels[temp], scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(
    x = "Chromosome",
    y = expression(-log[10](p-value)),
    title = "CMH Selection Scan Manhattan Plot (FDR < 0.05, Capped at 60)"
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("CMH_Manhattan_combined_qval_clipped.png", q, width = 12, height = 8, dpi = 300)





####plot all
q = ggplot(all_data, aes(x = global_pos, y = logp)) +
  geom_point(aes(color = CHROM), size = 0.5) +
  scale_color_manual(values = chrom_colors, guide = "none") +
  # Highlight significant points (qval < 0.05) in green
  geom_point(data = subset(all_data, qval < 0.05), aes(x = global_pos, y = logp), color = "red", size = 0.8) +
  scale_x_continuous(
    breaks = chrom_lengths$offset + chrom_lengths$max_pos / 2,
    labels = chrom_lengths$CHROM
  ) +
  facet_grid(pollinator_labels[pollinator] ~ temperature_labels[temp], scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(
    x = "Chromosome",
    y = expression(-log[10](p-value)),
    title = "CMH Selection Scan Manhattan Plot"
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("CMH_Manhattan_combined_qval.png", q, width = 12, height = 8, dpi = 300)



