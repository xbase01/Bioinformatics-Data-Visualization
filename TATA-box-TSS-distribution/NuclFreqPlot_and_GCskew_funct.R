library(Biostrings)
library(ggplot2)
library(seqinr)
library(stringr)
library(scales)
library(gridExtra)

setwd("C:/Users/judo4/Documents/Data_Visualization")

# ==============================================================================
# PART 1: NUCLEOTIDE FREQUENCY ANALYSIS
# ==============================================================================

file <- 'Pita.TSS_centered.1200.fasta' 
fasta <- readDNAStringSet(file, format="fasta")

# Calculate nucleotide frequencies
afmc <- consensusMatrix(fasta, as.prob=TRUE, baseOnly=TRUE)
tafmc <- as.data.frame(t(afmc))

# Add TSS-relative positions (-1000 to +199)
tafmc$pos <- seq(-1000, 199)
row.names(tafmc) <- tafmc$pos
tafmc$pos <- NULL

# Reshape for ggplot
rtafmc <- reshape2::melt(as.matrix(tafmc), value.name = 'Probability')
colnames(rtafmc) <- c('Position', 'Base', 'Probability')

# Define color scheme for bases (standard)
base_colors <- c("A" = "#5050FF", "C" = "#E6AB02", 
                 "G" = "#E7298A", "T" = "#66A61E")

# ==============================================================================
# PLOT 1: Full nucleotide composition across entire region
# ==============================================================================

p1 <- ggplot(rtafmc, aes(Position, Probability, color = Base)) +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E74C3C", 
             linewidth = 0.6, alpha = 0.7) +
  annotate("text", x = 20, y = max(rtafmc$Probability) * 0.95, 
           label = "TSS", color = "#E74C3C", size = 3.5, 
           fontface = "bold", hjust = 0) +
  scale_color_manual(values = base_colors,
                     labels = c("Adenine", "Cytosine", "Guanine", "Thymine")) +
  scale_x_continuous(breaks = seq(-1000, 200, 200),
                     labels = comma) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    x = "Position relative to TSS (bp)",
    y = "Nucleotide frequency",
    title = "Nucleotide Composition Across Promoter Region",
    color = "Base"
  ) +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5,
                              margin = margin(b = 15)),
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "black", 
                                      linewidth = 0.3),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 15, 10, 10)
  )

print(p1)

# Save full plot
ggsave("nucleotide_composition_full.png", plot = p1, 
       width = 8, height = 5, dpi = 600, bg = "white")
ggsave("nucleotide_composition_full.pdf", plot = p1, 
       width = 8, height = 5, device = cairo_pdf)
ggsave("nucleotide_composition_full.tiff", plot = p1, 
       width = 8, height = 5, dpi = 600, compression = "lzw")

# ==============================================================================
# PLOT 2: Zoomed view of promoter region (-80 to +20 bp)
# ==============================================================================

rtafmc_zoom <- subset(rtafmc, Position > -80 & Position < 20)

p2 <- ggplot(rtafmc_zoom, aes(Position, Probability, color = Base)) +
  geom_line(linewidth = 1.2, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E74C3C", 
             linewidth = 0.7, alpha = 0.8) +
  annotate("text", x = 2, y = max(rtafmc_zoom$Probability) * 0.95, 
           label = "TSS", color = "#E74C3C", size = 4, 
           fontface = "bold", hjust = 0) +
  # Highlight TATA box region
  annotate("rect", xmin = -35, xmax = -25, 
           ymin = 0, ymax = Inf,
           alpha = 0.1, fill = "gray") +
  annotate("text", x = -30, y = max(rtafmc_zoom$Probability) * 0.1, 
           label = "TATA\nregion", color = "gray30", size = 2.5, 
           hjust = 0.5, lineheight = 0.9) +
  scale_color_manual(values = base_colors,
                     labels = c("Adenine", "Cytosine", "Guanine", "Thymine")) +
  scale_x_continuous(breaks = seq(-80, 20, 20),
                     labels = comma) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    x = "Position relative to TSS (bp)",
    y = "Nucleotide frequency",
    title = "Nucleotide Composition in Core Promoter Region",
    subtitle = "Region: -80 to +20 bp relative to TSS",
    color = "Base"
  ) +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5,
                              margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, hjust = 0.5, 
                                  margin = margin(b = 10), face = "italic"),
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = "black", 
                                      linewidth = 0.3),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 15, 10, 10)
  )

print(p2)

# Save zoomed plot
ggsave("nucleotide_composition_zoom.png", plot = p2, 
       width = 8, height = 5, dpi = 600, bg = "white")
ggsave("nucleotide_composition_zoom.pdf", plot = p2, 
       width = 8, height = 5, device = cairo_pdf)
ggsave("nucleotide_composition_zoom.tiff", plot = p2, 
       width = 8, height = 5, dpi = 600, compression = "lzw")

# ==============================================================================
# PART 2: GC-SKEW ANALYSIS
# ==============================================================================

# Custom function to split sequences into windows
split.windows <- function(sequences, window_size = 40, step_size = 10) {
  n_seq <- length(sequences)
  seq_length <- nchar(sequences[1])
  
  # Calculate number of windows
  n_windows <- floor((seq_length - window_size) / step_size) + 1
  
  # Initialize matrix
  window_matrix <- matrix(NA, nrow = n_seq, ncol = n_windows)
  
  for (i in 1:n_windows) {
    start <- 1 + (i - 1) * step_size
    end <- start + window_size - 1
    
    for (j in 1:n_seq) {
      window_matrix[j, i] <- substring(as.character(sequences[j]), start, end)
    }
  }
  
  return(data.frame(window_matrix))
}

# Function to compute GC-skew
GCskew <- function(fasta, window_size = 40, step_size = 10) {
  
  windows <- split.windows(as.vector(fasta), window_size, step_size)
  
  # Count C and G in each window
  c_counts <- apply(windows, MARGIN = 2, FUN = function(x) {
    sum(str_count(x, 'C'))
  })
  
  g_counts <- apply(windows, MARGIN = 2, FUN = function(x) {
    sum(str_count(x, 'G'))
  })
  
  # Compute GC-skew: (G-C)/(G+C)
  skew <- (g_counts - c_counts) / (g_counts + c_counts)
  
  # Calculate window start positions
  seq_length <- nchar(fasta[1])
  n_windows <- floor((seq_length - window_size) / step_size) + 1
  window_starts <- seq(1, 1 + (n_windows - 1) * step_size, by = step_size)
  
  # Convert to TSS-relative coordinates
  # Position 1 in sequence = -1000 bp relative to TSS
  relative_positions <- window_starts - 1001
  
  skew_df <- data.frame(
    Position = relative_positions,
    GCskew = skew
  )
  
  return(skew_df)
}

# Apply the function
w <- 40   # window length
s <- 10   # step size

cat("\n=== Computing GC-skew ===\n")
cat("Window size:", w, "bp\n")
cat("Step size:", s, "bp\n")

df_skew <- GCskew(fasta, window_size = w, step_size = s)

cat("Total windows analyzed:", nrow(df_skew), "\n")
cat("Position range:", min(df_skew$Position), "to", max(df_skew$Position), "bp\n\n")

# ==============================================================================
# PLOT 3: GC-skew across entire region
# ==============================================================================

p3 <- ggplot(df_skew, aes(Position, GCskew)) +
  geom_line(linewidth = 0.8, color = "#2C3E50") +
  geom_hline(yintercept = 0, linetype = "solid", 
             color = "gray40", linewidth = 0.4, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E74C3C", 
             linewidth = 0.6, alpha = 0.7) +
  annotate("text", x = 20, y = max(df_skew$GCskew) * 0.9, 
           label = "TSS", color = "#E74C3C", size = 3.5, 
           fontface = "bold", hjust = 0) +
  # Add labels for positive/negative skew
  annotate("text", x = -900, y = max(df_skew$GCskew) * 0.9, 
           label = "G-rich", color = "gray30", size = 3, 
           fontface = "italic", hjust = 0) +
  annotate("text", x = -900, y = min(df_skew$GCskew) * 0.9, 
           label = "C-rich", color = "gray30", size = 3, 
           fontface = "italic", hjust = 0) +
  scale_x_continuous(breaks = seq(-1000, 200, 200),
                     labels = comma) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    x = "Position relative to TSS (bp)",
    y = "GC-skew [(G-C)/(G+C)]",
    title = "GC-Skew Analysis Across Promoter Region",
    subtitle = paste0("Window size: ", w, " bp, Step size: ", s, " bp")
  ) +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5,
                              margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, hjust = 0.5, 
                                  margin = margin(b = 10), face = "italic"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 15, 10, 10)
  )

print(p3)

# Save GC-skew plot
ggsave("GC_skew_analysis.png", plot = p3, 
       width = 8, height = 5, dpi = 600, bg = "white")
ggsave("GC_skew_analysis.pdf", plot = p3, 
       width = 8, height = 5, device = cairo_pdf)
ggsave("GC_skew_analysis.tiff", plot = p3, 
       width = 8, height = 5, dpi = 600, compression = "lzw")

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

cat("=== Nucleotide Composition Summary ===\n")
cat("Total sequences analyzed:", length(fasta), "\n")
cat("Sequence length:", nchar(fasta[1]), "bp\n")
cat("Coordinate range: -1000 to +199 bp (TSS at position 0)\n\n")

# Calculate average nucleotide frequencies
avg_freqs <- colMeans(tafmc)
cat("Average nucleotide frequencies:\n")
cat(sprintf("  A: %.2f%%\n", avg_freqs["A"] * 100))
cat(sprintf("  C: %.2f%%\n", avg_freqs["C"] * 100))
cat(sprintf("  G: %.2f%%\n", avg_freqs["G"] * 100))
cat(sprintf("  T: %.2f%%\n", avg_freqs["T"] * 100))
cat(sprintf("  GC content: %.2f%%\n", (avg_freqs["G"] + avg_freqs["C"]) * 100))

cat("\n=== GC-Skew Summary ===\n")
cat("Mean GC-skew:", round(mean(df_skew$GCskew), 4), "\n")
cat("Median GC-skew:", round(median(df_skew$GCskew), 4), "\n")
cat("Min GC-skew:", round(min(df_skew$GCskew), 4), 
    " at position", df_skew$Position[which.min(df_skew$GCskew)], "bp\n")
cat("Max GC-skew:", round(max(df_skew$GCskew), 4), 
    " at position", df_skew$Position[which.max(df_skew$GCskew)], "bp\n")

cat("\n=== Output Files Generated ===\n")
cat("Nucleotide composition (full): nucleotide_composition_full.[png/pdf/tiff]\n")
cat("Nucleotide composition (zoom): nucleotide_composition_zoom.[png/pdf/tiff]\n")
cat("GC-skew analysis: GC_skew_analysis.[png/pdf/tiff]\n")
cat("\nAll figures saved at 600 DPI in multiple formats.\n")