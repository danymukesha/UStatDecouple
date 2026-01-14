# Analyze Gene Expression Correlation Structure

Case study analyzing correlation structure in gene expression data

## Usage

``` r
analyze_gene_expression_correlations(
  num_genes = 20,
  num_samples = 15,
  B = 500,
  seed = 123
)
```

## Arguments

- num_genes:

  Number of genes to simulate

- num_samples:

  Number of samples/conditions

- B:

  Number of decoupling iterations

- seed:

  Random seed for reproducibility

## Value

DecoupleResult object with analysis results
