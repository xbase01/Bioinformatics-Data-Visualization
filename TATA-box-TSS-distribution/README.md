# TATA Box Distribution Analysis in Rice Promoters

[![R](https://img.shields.io/badge/R-4.0%2B-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This project analyzes where TATA boxes appear in rice gene promoters. TATA boxes are conserved DNA sequences found in the core promoter region of many eukaryotic genes. They act as a recognition site for the TATA-binding protein (TBP), which is a component of the general transcription factor TFIID.

## What We're Studying

**Research Question**: Where exactly do TATA boxes sit relative to where transcription starts (the TSS)?

**Why It Matters**: TATA boxes need to be in the right place to work properly. Understanding their positioning helps us predict gene behavior and engineer better crops.

## The Data

- **Source**: Rice (*Oryza sativa*) gene promoter sequences
- **Dataset**: 9,180 promoter sequences
- **Structure**: Each sequence is 1,200 bp long, centered on the transcription start site
  - Upstream region: -1000 to 0 bp
  - Downstream region: 0 to +199 bp
- **Strands**: Includes both plus (+) and minus (-) strand genes

## Quick Start

### 1. Install Required Packages

```r
# Install Bioconductor manager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install packages
BiocManager::install("Biostrings")
install.packages(c("ggplot2", "scales"))
```

### 2. Run the Analysis

```r
# Set your working directory
setwd("path/to/your/folder")

# Run the script
source("tata_box_analysis.R")
```

### 3. Check Your Results

The script generates three publication-ready figures:
- `TATA_distribution_TSS.png` (600 DPI)
- `TATA_distribution_TSS.pdf` (vector format)
- `TATA_distribution_TSS.tiff` (journal submission)

## How It Works

1. **Load sequences**: Reads all 9,180 rice promoter sequences
2. **Fix strand orientation**: Reverse complements minus-strand sequences so all read in the same direction
3. **Sliding window search**: Uses a 40 bp window that slides in 10 bp steps, looking for "TATA" sequences
4. **Calculate frequency**: Determines what percentage of genes have TATA at each position
5. **Plot results**: Creates a graph showing TATA box distribution

## Key Parameters

```r
motif <- "TATA"           # What we're searching for
window_width <- 40        # Search window size (bp)
window_increment <- 10    # Step size (bp)
tss_position <- 1001      # Where TSS is in the sequence
```

## Understanding Your Results

### What to Look For

✅ **Peak around -20 to -35 bp**: This is where TATA boxes should be  
✅ **~40% frequency**: Typical for rice promoters  
✅ **Sharp peak upstream**: Shows TATA boxes are precisely positioned  

### Example Results

From our analysis of 9,180 rice genes:
- **Peak position**: -21 bp (right where expected!)
- **Peak frequency**: 40% of genes
- **Canonical region (-35 to -25 bp)**: 39.85% average frequency
- **Estimated functional TATA boxes**: ~3,658 genes

## File Structure

```
TATA-box-TSS-distribution/
├── README.md                                  # This file
├── tata_box_analysis.R                        # Main analysis script
├── Pita.TSS_centered.1200.fasta              # Rice promoter sequences
└── results/
    ├── TATA_distribution_TSS.png             # Figure (600 DPI)
    ├── TATA_distribution_TSS.pdf             # Figure (vector)
    └── TATA_distribution_TSS.tiff            # Figure (TIFF)
    ...
```

## Customizing the Analysis

### Search for Different Motifs

```r
motif <- "CAAT"   # CAAT box
motif <- "CCAAT"  # Extended CAAT box
motif <- "GC"     # GC box
```

### Adjust Window Size

```r
window_width <- 50        # Larger window (smoother curve)
window_increment <- 5     # Smaller steps (finer resolution)
```

## Common Issues

**Q: Plot looks empty or wrong**  
A: Check that your FASTA file is in the working directory and properly formatted

**Q: Peak is in the wrong place**  
A: Make sure `tss_position` matches your sequence structure (1001 for these sequences)

**Q: Package errors**  
A: Install Bioconductor packages first, then CRAN packages

## What's Special About This Analysis

- ✅ **Strand-corrected**: Properly handles genes on both DNA strands
- ✅ **Publication-ready**: Generates high-quality figures (600 DPI)
- ✅ **Biologically validated**: Results match expected TATA box positioning
- ✅ **Large-scale**: Analyzes thousands of sequences in minutes

## Beyond TATA Boxes

This same approach works for:
- 🔍 Other promoter elements (CAAT box, GC box, Inr)
- 🔍 Response elements (stress, hormone, light)
- 🔍 Splice sites and polyadenylation signals
- 🔍 Any DNA motif you want to map


## License

MIT License - Free to use for research and education

## Questions?

Open an issue on GitHub or contact the repository maintainer.

---

**Built with**: R, Biostrings, ggplot2  
**Last Updated**: December 2025  
**Status**: Active