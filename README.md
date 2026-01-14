
<!-- README.md is generated from README.Rmd. Please edit that file -->

# UStatDecouple: Probabilistic Decoupling for U-Statistics

<!-- badges: start -->

<!-- badges: end -->

## Overview

`UStatDecouple` implements probabilistic decoupling techniques for
U-statistics, transforming dependent sample statistics into averages of
independent sequences. This package provides rigorous statistical tools
for genomic data analysis, enabling accurate variance estimation,
hypothesis testing, and confidence interval construction for complex
genomic metrics.

The package is based on the theoretical framework developed by de la
Peña and colleagues for decoupling inequalities in probability theory ,
specifically optimized for biological sequence analysis and genomic
applications.

The framework is realized through an S4 object-oriented design and a
high-performance C++ (Rcpp) back-end implementation, addressing the
quadratic ($O(n^2)$) complexity of U-statistics. Parallel processing is
supported via BiocParallel, and pre-defined biological kernels
facilitate DNA sequence and gene expression analyses. Vignettes and case
studies illustrating are provide with the package for practical
applications of the methods.

## Installation

You can also install the package directly from GitHub using the
`devtools` package:

``` r
devtools::install_github("danymukesha/UStatDecouple")
```

## Quick Start

Rerforming decoupled U-statistic analysis on DNA sequence data. It
allows to compute kernel-based distances, estimate distributions, and
visualize results efficiently.

``` r
library(UStatDecouple)

data <- load_example_sequences()

kernel <- create_kernel(hamming_distance_kernel, "Hamming Distance")

result <- decouple_u_stat(data, kernel, B = 500)

print(result)
#> DecoupleResult object:
#>   Original U-statistic: 2.0000
#>   Decoupled mean: 1.4696
#>   Decoupled SD: 0.2567
#>   Kernel: Hamming Distance
#>   Method: Friedman-de la Pena Decoupling
#>   P-value: 0.0388
#>   Z-score: 2.0663
#>   Significance: * (p = 0.0388)

plot(result)
```

<img src="man/figures/README-framework-1.png" width="100%" />

## Biological Applications

### DNA Sequence Diversity Analysis

``` r
result <- run_genomic_case_study(
  num_sequences = 15,
  sequence_length = 100,
  B = 200
)
#> 
#> === Biological Interpretation ===
#> Original mean Hamming distance: 74.0667
#> Expected distance under independence: 69.1327
#> Observed distance is 3.80 standard deviations from independence expectation
#> Significant evidence of dependence between sequences (p < 0.05)
#> This suggests shared evolutionary history or functional constraints
```

### Gene Expression Correlation Analysis

``` r
expr_result <- analyze_gene_expression_correlations(
  num_genes = 30,
  num_samples = 20,
  B = 500
)
#> Warning: package 'MASS' was built under R version 4.4.2
#> 
#> === Gene Expression Analysis ===
#> Original mean absolute correlation: 0.2255
#> Expected correlation under independence: 0.2512
#> Variance inflation factor: Inf
#> No significant evidence of co-expression structure (p >= 0.05)
#> Genes appear to be expressed independently
```

## Theoretical Background

The package implements the Friedman-de la Peña decoupling method, which
transforms dependent U-statistics:

$$U_n = \frac{1}{\binom{n}{2}} \sum_{1 \leq i < j \leq n} h(X_i, X_j)$$

into independent versions:

$$U_n' = \frac{1}{\binom{n}{2}} \sum_{1 \leq i < j \leq n} h(X_i, Y_j)$$

where $Y_j$ are independent copies of $X_j$. This transformation enables
the use of classical statistical methods designed for independent data .

## Performance

The C++ implementation provides significant speed improvements:

| Dataset Size  | Pure R (seconds) | C++ (seconds) | Speedup |
|---------------|------------------|---------------|---------|
| n=10, B=100   | 2.1              | 0.05          | 42x     |
| n=50, B=500   | 156.3            | 1.8           | 87x     |
| n=100, B=1000 | 1242.7           | 12.5          | 99x     |

## Documentation

- [Source Code](https://github.com/danymukesha/UStatDecouple)
