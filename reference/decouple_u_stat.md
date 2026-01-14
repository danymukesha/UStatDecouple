# Decouple a U-statistic

Implements probabilistic decoupling for U-statistics using independent
sequence copies

## Usage

``` r
decouple_u_stat(x, kernel, B = 1000, parallel = FALSE, seed = 123)
```

## Arguments

- x:

  List of data points (sequences, expression profiles, etc.)

- kernel:

  UStatKernel object or function

- B:

  Number of decoupling iterations

- parallel:

  Logical indicating whether to use parallel computation

- seed:

  Random seed for reproducibility

## Value

DecoupleResult object containing original and decoupled statistics

## Examples

``` r
# Generate synthetic DNA sequences
set.seed(123)
bases <- c("A", "C", "G", "T")
sequences <- lapply(1:5, function(i) sample(bases, 20, replace = TRUE))

# Create kernel and run decoupling
kernel <- create_kernel(hamming_distance_kernel, "Hamming Distance")
result <- decouple_u_stat(sequences, kernel, B = 100)

# Print results
result
#> DecoupleResult object:
#>   Original U-statistic: 15.4000
#>   Decoupled mean: 12.0850
#>   Decoupled SD: 1.8846
#>   Kernel: Hamming Distance
#>   Method: Friedman-de la Pena Decoupling
#>   P-value: 0.0786
#>   Z-score: 1.7590
#>   Significance:   (p = 0.0786)
```
