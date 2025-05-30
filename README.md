# TAD Consensus Analysis

## Overview

This R package provides tools for generating consensus Topologically Associating Domains (TADs) from multiple prediction
methods. TADs are fundamental units of chromatin organization that play crucial roles in gene regulation. While multiple
computational tools exist to predict TAD boundaries from Hi-C data, their results often vary significantly. This package
implements methods to integrate predictions from multiple tools and generate high-confidence consensus TAD sets.

## Installation

```r
# Install from GitHub
devtools::install_github("CSOgroup/consensusTADs", build_vignettes = TRUE)
```

## Key Features

- Generate consensus TADs from multiple prediction tools
- Calculate Measure of Concordance (MoC) between TAD predictions
- Select optimal non-overlapping TAD sets using dynamic programming
- Apply iterative threshold approach for consensus building

## Main Functions

### `generate_tad_consensus()`

Creates consensus TADs through an iterative threshold approach that selects optimal non-overlapping TADs representing
agreement across different prediction methods.

```r
consensus_tads <- generate_tad_consensus(
  df_tools,      # Data frame with TAD predictions
  threshold = 0, # Minimum MoC threshold
  step = -0.05   # Step size for threshold iteration
)
```

### `moc_score_filter()`

Calculates the Measure of Concordance (MoC) between TAD predictions and filters significant overlaps based on a
threshold.

### `select_global_optimal_tads()`

Implements a dynamic programming algorithm to select a set of non-overlapping TADs that maximize the total MoC score.

## Example Usage

```r
# Prepare input data with predictions from multiple tools
tad_data <- data.frame(
  chr = rep("chr1", 6),
  start = c(10000, 20000, 50000, 12000, 22000, 48000),
  end = c(30000, 45000, 65000, 32000, 43000, 67000),
  meta.tool = c(rep("tool1", 3), rep("tool2", 3))
)

# Generate consensus TADs with default parameters
library(consensusTADs)
consensus_results <- generate_tad_consensus(tad_data)
print(consensus_results)

# Generate consensus TADs with custom threshold values
custom_consensus <- generate_tad_consensus(
  tad_data,
  threshold = 0.3,
  step = -0.1
)
```

## How It Works

The consensus generation process follows these steps:

1. **Input validation**: Check if the input contains data from multiple prediction tools
2. **Data preparation**: Split the input data by chromosome
3. **Threshold sequence generation**: Create a sequence of threshold values
4. **Iterative TAD selection**: For each chromosome and threshold, calculate MoC scores and select optimal TADs
5. **Result compilation**: Combine results from all chromosomes

## The Measure of Concordance (MoC) Score

The MoC score quantifies the agreement between two TAD predictions:

MoC = (intersection_width)² / (width1 × width2)

Where:

- `intersection_width` is the length of the overlap between two TADs
- `width1` and `width2` are the lengths of the two TADs being compared

## Dependencies

- dplyr
- GenomicRanges
- IRanges
- tibble
- purrr
- tidyr
- stringr
- magrittr

