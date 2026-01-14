#' @title Run Genomic Case Study: DNA Sequence Diversity
#' @description Case study analyzing genetic diversity using decoupling of Hamming distances
#' @param num_sequences Number of DNA sequences to generate
#' @param sequence_length Length of each DNA sequence
#' @param B Number of decoupling iterations
#' @param seed Random seed for reproducibility
#' @return DecoupleResult object with analysis results
#' @export
#' @examples
#' # Run case study with default parameters
#' result <- run_genomic_case_study()
#'
#' # Plot results
#' plot(result)
run_genomic_case_study <- function(num_sequences = 10, sequence_length = 50,
                                   B = 500, seed = 123) {
    set.seed(seed)

    # Generate synthetic but biologically realistic DNA sequences
    bases <- c("A", "C", "G", "T")

    # Create sequences with some biological structure (not completely random)
    sequences <- vector("list", num_sequences)

    for (i in 1:num_sequences) {
        # Start with a random sequence
        seq <- sample(bases, sequence_length, replace = TRUE, prob = c(0.3, 0.2, 0.3, 0.2))

        # Introduce some correlation structure (simulating evolutionary relationships)
        if (i > 1) {
            mutation_rate <- 0.1 # 10% mutation rate from previous sequence
            mutation_positions <- sample(sequence_length, size = round(sequence_length * mutation_rate))
            seq[mutation_positions] <- sample(bases, size = length(mutation_positions), replace = TRUE)
        }

        sequences[[i]] <- seq
    }

    # Create Hamming distance kernel
    hamming_kernel <- create_kernel(hamming_distance_kernel, "Hamming Distance")

    # Run decoupling analysis
    result <- decouple_u_stat(sequences, hamming_kernel, B = B, seed = seed)

    # Add biological interpretation
    cat("\n=== Biological Interpretation ===\n")
    cat(sprintf("Original mean Hamming distance: %.4f\n", result@original_stat))
    cat(sprintf("Expected distance under independence: %.4f\n", mean(result@decoupled_distribution)))
    cat(sprintf(
        "Observed distance is %.2f standard deviations from independence expectation\n",
        (result@original_stat - mean(result@decoupled_distribution)) / sd(result@decoupled_distribution)
    ))

    if (result@p_value < 0.05) {
        cat("Significant evidence of dependence between sequences (p < 0.05)\n")
        cat("This suggests shared evolutionary history or functional constraints\n")
    } else {
        cat("No significant evidence of dependence between sequences (p >= 0.05)\n")
        cat("Sequences appear to evolve independently\n")
    }

    return(result)
}

#' @title Analyze Gene Expression Correlation Structure
#' @description Case study analyzing correlation structure in gene expression data
#' @param num_genes Number of genes to simulate
#' @param num_samples Number of samples/conditions
#' @param B Number of decoupling iterations
#' @param seed Random seed for reproducibility
#' @return DecoupleResult object with analysis results
#' @export
analyze_gene_expression_correlations <- function(num_genes = 20, num_samples = 15,
                                                 B = 500, seed = 123) {
    set.seed(seed)

    # Simulate gene expression data with biological structure
    # Some genes will be co-expressed (dependent), others independent

    # Create correlation matrix with block structure (simulating pathways)
    block_size <- 5
    num_blocks <- num_genes %/% block_size

    # Initialize correlation matrix
    cor_matrix <- diag(num_genes)

    # Add correlation within blocks
    for (block in 1:num_blocks) {
        start_idx <- (block - 1) * block_size + 1
        end_idx <- min(block * block_size, num_genes)

        # High correlation within block (simulating co-regulated genes)
        for (i in start_idx:end_idx) {
            for (j in start_idx:end_idx) {
                if (i != j) {
                    cor_matrix[i, j] <- 0.7 + runif(1, -0.1, 0.1) # Correlation around 0.7
                }
            }
        }
    }

    # Make correlation matrix positive definite
    eigen_vals <- eigen(cor_matrix)$values
    if (min(Re(eigen_vals)) < 0) {
        cor_matrix <- cor_matrix + diag(abs(min(eigen_vals)) + 0.1, num_genes)
    }

    # Generate multivariate normal data
    library(MASS)
    expr_data <- mvrnorm(n = num_samples, mu = rep(0, num_genes), Sigma = cor_matrix)
    expr_data <- t(data.frame(expr_data))

    # Convert to expression profiles (rows = genes, columns = samples)
    rownames(expr_data) <- paste0("Gene_", 1:num_genes)
    colnames(expr_data) <- paste0("Sample_", 1:num_samples)

    # Create list of gene expression profiles
    gene_profiles <- split(expr_data, row(expr_data))

    # Create correlation kernel
    corr_kernel <- create_kernel(
        gene_expression_correlation_kernel,
        "Absolute Spearman Correlation"
    )

    # Run decoupling analysis
    result <- decouple_u_stat(gene_profiles, corr_kernel, B = B, seed = seed)

    # Add biological interpretation
    cat("\n=== Gene Expression Analysis ===\n")
    cat(sprintf("Original mean absolute correlation: %.4f\n", result@original_stat))
    cat(sprintf("Expected correlation under independence: %.4f\n", mean(result@decoupled_distribution)))
    cat(sprintf(
        "Variance inflation factor: %.2f\n",
        var(result@decoupled_distribution) / var(rep(result@original_stat, B))
    ))

    if (result@p_value < 0.05) {
        cat("Significant evidence of co-expression structure (p < 0.05)\n")
        cat("This suggests the presence of regulatory networks or functional modules\n")
    } else {
        cat("No significant evidence of co-expression structure (p >= 0.05)\n")
        cat("Genes appear to be expressed independently\n")
    }

    return(result)
}

#' @title Plot DecoupleResult
#' @param x DecoupleResult object
#' @param main Plot title
#' @param ... Additional arguments to plot
#' @export
setMethod("plot", "DecoupleResult", function(x, main = "Decoupling Analysis", ...) {
    # Create histogram of decoupled distribution
    hist(x@decoupled_distribution,
        breaks = 30,
        col = "lightblue",
        border = "white",
        main = main,
        xlab = "Decoupled Statistic Value",
        ...
    )

    # Add vertical line for original statistic
    abline(v = x@original_stat, col = "red", lwd = 2, lty = 2)
    legend("topright",
        legend = sprintf("Original: %.4f", x@original_stat),
        col = "red", lwd = 2, lty = 2, bty = "n"
    )

    # Add text for p-value
    text(0.95 * par("usr")[2], 0.95 * par("usr")[4],
        sprintf("p = %.4f", x@p_value),
        pos = 2, cex = 0.8
    )
})
