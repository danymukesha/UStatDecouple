# Run Genomic Case Study: DNA Sequence Diversity

Case study analyzing genetic diversity using decoupling of Hamming
distances

## Usage

``` r
run_genomic_case_study(
  num_sequences = 10,
  sequence_length = 50,
  B = 500,
  seed = 123
)
```

## Arguments

- num_sequences:

  Number of DNA sequences to generate

- sequence_length:

  Length of each DNA sequence

- B:

  Number of decoupling iterations

- seed:

  Random seed for reproducibility

## Value

DecoupleResult object with analysis results

## Examples

``` r
# Run case study with default parameters
result <- run_genomic_case_study()
#> 
#> === Biological Interpretation ===
#> Original mean Hamming distance: 37.4889
#> Expected distance under independence: 33.6505
#> Observed distance is 3.46 standard deviations from independence expectation
#> Significant evidence of dependence between sequences (p < 0.05)
#> This suggests shared evolutionary history or functional constraints

# Plot results
plot(result)
```
