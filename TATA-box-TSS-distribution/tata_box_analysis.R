library(Biostrings)
library(ggplot2)
library(scales)

setwd("C:/Users/judo4/Documents/Data_Visualization")

# ==============================================================================
# TATA BOX DISTRIBUTION ANALYSIS
# Considered for biological relevance: (1) TSS position, (2) Strand orientation, (3) Consistent windows
# ==============================================================================

# Read multi-fasta file
sequences <- readDNAStringSet("Pita.TSS_centered.1200.fasta", format = "fasta")

cat("=== Processing Strand Orientations ===\n")
# Extract strand information from sequence names
seq_names <- names(sequences)
is_minus_strand <- grepl("\\(-\\)", seq_names)

cat("Plus strand sequences:", sum(!is_minus_strand), "\n")
cat("Minus strand sequences:", sum(is_minus_strand), "\n\n")

# Reverse complement minus strand sequences so all are in same orientation
for (i in which(is_minus_strand)) {
  sequences[[i]] <- reverseComplement(sequences[[i]])
}

cat("All sequences now oriented 5' to 3' relative to TSS\n\n")

# Function to count motif occurrences in a sliding window
CountMotifSlidingWindow <- function(sequence, motif, window_width, window_increment) {
  seq_length <- nchar(sequence)
  n_windows <- (seq_length - window_width) %/% window_increment + 1
  motif_counts <- integer(n_windows)
  
  for (i in seq_along(motif_counts)) {
    start <- 1 + (i - 1) * window_increment
    end <- start + window_width - 1
    window <- substring(as.character(sequence), start, end)
    motif_counts[i] <- sum(gregexpr(motif, window, fixed = TRUE)[[1]] > 0)
  }
  
  return(motif_counts)
}

# Set parameters - Consistent with NuclFreqPlot and GC skew analysis
motif <- "TATA"
window_width <- 40    # Matching GC-skew analysis
window_increment <- 10 # Matching other analyses
tss_position <- 1001   # TSS is at position 1001 (maps to 0 in plot)

cat("=== Analysis Parameters ===\n")
cat("Motif searched:", motif, "\n")
cat("Window width:", window_width, "bp (matching GC-skew analysis)\n")
cat("Window increment:", window_increment, "bp\n")
cat("TSS position in sequence:", tss_position, "(position 0 in plot)\n\n")

# Initialize matrix to store motif counts
seq_length <- nchar(sequences[[1]])
n_windows <- (seq_length - window_width) %/% window_increment + 1
motif_matrix <- matrix(NA, nrow = length(sequences), ncol = n_windows)

# Fill the matrix with motif counts for each sequence
cat("Analyzing", length(sequences), "sequences...\n")
for (i in seq_along(sequences)) {
  motif_matrix[i, ] <- CountMotifSlidingWindow(sequences[[i]], motif, 
                                                 window_width, window_increment)
  
  if (i %% 1000 == 0) {
    cat("  Processed", i, "sequences...\n")
  }
}

# Calculate the frequency of sequences with TATA motif in each window
window_sequence_counts <- colSums(motif_matrix > 0, na.rm = TRUE)
window_frequency <- (window_sequence_counts / nrow(motif_matrix)) * 100

# Calculate window centers relative to TSS
window_centers <- seq(window_width/2, seq_length - window_width/2, by = window_increment)
relative_positions <- window_centers - tss_position

# Create data frame for plotting
df <- data.frame(
  Position = relative_positions,
  Frequency = window_frequency,
  Count = window_sequence_counts
)

# Create quality plot
p <- ggplot(df, aes(x = Position, y = Frequency)) +
  geom_line(linewidth = 0.8, color = "#2C3E50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E74C3C", 
             linewidth = 0.6, alpha = 0.7) +
  annotate("text", x = 10, y = max(df$Frequency) * 0.95, label = "TSS", 
           color = "#E74C3C", size = 3.5, fontface = "bold", hjust = 0) +
  # Add shaded region for canonical TATA box location
  annotate("rect", xmin = -35, xmax = -25, ymin = 0, ymax = Inf,
           alpha = 0.15, fill = "blue") +
  annotate("text", x = -30, y = max(df$Frequency) * 0.05, 
           label = "Canonical", 
           color = "blue", size = 2.8, hjust = 0.5) +
  scale_x_continuous(breaks = seq(-1000, 200, 200),
                     labels = comma) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    x = "Position relative to TSS (bp)",
    y = "Frequency (%)",
    title = "TATA Box Distribution Relative to Transcription Start Site",
    subtitle = paste0("Window: ", window_width, " bp, Step: ", window_increment, 
                     " bp")
  ) +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    # Axis styling
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    
    # Title styling
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5,
                              margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, hjust = 0.5, 
                                  margin = margin(b = 10), face = "italic"),
    
    # Panel styling
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Plot margins
    plot.margin = margin(15, 15, 10, 10)
  )

# Display plot
print(p)

# Save high-resolution figures
ggsave("TATA_distribution_TSS.png", plot = p, 
       width = 8, height = 5, dpi = 600, bg = "white")

ggsave("TATA_distribution_TSS.pdf", plot = p, 
       width = 8, height = 5, device = cairo_pdf)

ggsave("TATA_distribution_TSS.tiff", plot = p, 
       width = 8, height = 5, dpi = 600, compression = "lzw")

# ==============================================================================
# DETAILED SUMMARY STATISTICS
# ==============================================================================

cat("\n=== Analysis Summary ===\n")
cat("Total sequences analyzed:", nrow(motif_matrix), "\n")
cat("  Plus strand:", sum(!is_minus_strand), "\n")
cat("  Minus strand (reverse complemented):", sum(is_minus_strand), "\n\n")

cat("Motif searched:", motif, "\n")
cat("Window width:", window_width, "bp\n")
cat("Window increment:", window_increment, "bp\n")
cat("Coordinate range: -980 to +179 bp relative to TSS\n\n")

cat("Peak frequency:", round(max(window_frequency), 2), "%\n")
cat("Peak position relative to TSS:", 
    df$Position[which.max(df$Frequency)], "bp\n\n")

# Additional biological context
peak_pos <- df$Position[which.max(df$Frequency)]
cat("=== Biological Interpretation ===\n")

if (peak_pos >= -50 && peak_pos <= -20) {
  cat("✓ Peak is in canonical TATA box region (-50 to -20 bp)\n")
  cat("  This suggests functional TATA boxes in rice promoters\n")
} else if (peak_pos < -20) {
  cat("! Peak is upstream of canonical position\n")
  cat("  Position:", peak_pos, "bp (canonical: -25 to -35 bp)\n")
  cat("  This may indicate:\n")
  cat("  - Broader distribution of functional TATA boxes in rice\n")
  cat("  - Alternative promoter architectures\n")
  cat("  - Presence of TATA-like sequences with regulatory functions\n")
} else {
  cat("! Peak is downstream of expected position\n")
  cat("  This requires further investigation\n")
}

# Calculate enrichment in different regions
upstream_strong <- df[df$Position >= -50 & df$Position <= -20, ]
upstream_extended <- df[df$Position >= -100 & df$Position <= -50, ]
downstream <- df[df$Position >= 0 & df$Position <= 100, ]

cat("\n=== Regional Analysis ===\n")
cat("Strong upstream region (-50 to -20 bp):\n")
cat("  Mean frequency:", round(mean(upstream_strong$Frequency), 2), "%\n")
cat("Extended upstream region (-100 to -50 bp):\n")
cat("  Mean frequency:", round(mean(upstream_extended$Frequency), 2), "%\n")
cat("Downstream region (0 to +100 bp):\n")
cat("  Mean frequency:", round(mean(downstream$Frequency), 2), "%\n")

# Find sequences with TATA in canonical region
canonical_window <- which(df$Position >= -35 & df$Position <= -25)
canonical_freq <- mean(df$Frequency[canonical_window])
cat("\nCanonical TATA region (-35 to -25 bp):\n")
cat("  Average frequency:", round(canonical_freq, 2), "%\n")
cat("  Estimated", round(canonical_freq * length(sequences) / 100), 
    "sequences have TATA in canonical position\n")

cat("\n=== Output Files ===\n")
cat("Figures saved:\n")
cat("  - TATA_distribution_TSS.png (600 DPI)\n")
cat("  - TATA_distribution_TSS.pdf (vector)\n")
cat("  - TATA_distribution_TSS.tiff (600 DPI, LZW compressed)\n")