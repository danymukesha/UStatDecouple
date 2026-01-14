#' @title Hamming Distance Kernel
#' @description Kernel function for calculating Hamming distance between DNA sequences
#' @param seq1 First DNA sequence (character vector)
#' @param seq2 Second DNA sequence (character vector)
#' @return Numeric Hamming distance
#' @export
hamming_distance_kernel <- function(seq1, seq2) {
    if (length(seq1) != length(seq2)) {
        stop("Sequences must be of equal length")
    }
    sum(seq1 != seq2)
}

#' @title Gene Expression Correlation Kernel
#' @description Kernel function for calculating absolute correlation between gene expression profiles
#' @param expr1 First gene expression profile (numeric vector)
#' @param expr2 Second gene expression profile (numeric vector)
#' @return Numeric absolute correlation
#' @export
gene_expression_correlation_kernel <- function(expr1, expr2) {
    if (length(expr1) != length(expr2)) {
        stop("Expression profiles must be of equal length")
    }
    abs(cor(expr1, expr2, method = "spearman"))
}

#' @title Decouple a U-statistic
#' @description Implements probabilistic decoupling for U-statistics using independent sequence copies
#' @param x List of data points (sequences, expression profiles, etc.)
#' @param kernel UStatKernel object or function
#' @param B Number of decoupling iterations
#' @param parallel Logical indicating whether to use parallel computation
#' @param seed Random seed for reproducibility
#' @return DecoupleResult object containing original and decoupled statistics
#' @importFrom stats sd
#' @importFrom BiocParallel bplapply MulticoreParam
#' @export
#' @examples
#' # Generate synthetic DNA sequences
#' set.seed(123)
#' bases <- c("A", "C", "G", "T")
#' sequences <- lapply(1:5, function(i) sample(bases, 20, replace = TRUE))
#'
#' # Create kernel and run decoupling
#' kernel <- create_kernel(hamming_distance_kernel, "Hamming Distance")
#' result <- decouple_u_stat(sequences, kernel, B = 100)
#'
#' # Print results
#' result
decouple_u_stat <- function(x, kernel, B = 1000, parallel = FALSE, seed = 123) {
    # Input validation
    if (length(x) < 2) stop("Sample size must be at least 2 for U-statistics")

    # Handle kernel input
    if (inherits(kernel, "UStatKernel")) {
        kernel_function <- kernel@kernel_function
        kernel_name <- kernel@kernel_name
    } else if (is.function(kernel)) {
        kernel_function <- kernel
        kernel_name <- "Custom Kernel"
    } else {
        stop("Kernel must be either a UStatKernel object or a function")
    }

    n <- length(x)
    set.seed(seed)

    # Calculate original U-statistic
    if (kernel@symmetric) {
        # For symmetric kernels, use combinations without replacement
        pairs <- combn(seq_len(n), 2)
        original_values <- apply(pairs, 2, function(idx) {
            kernel_function(x[[idx[1]]], x[[idx[2]]])
        })
        original_stat <- mean(original_values)
    } else {
        # For asymmetric kernels, use all pairs with i != j
        original_values <- numeric(n * (n - 1))
        count <- 1
        for (i in seq_len(n)) {
            for (j in seq_len(n)) {
                if (i != j) {
                    original_values[count] <- kernel_function(x[[i]], x[[j]])
                    count <- count + 1
                }
            }
        }
        original_stat <- mean(original_values)
    }

    # Setup parallel processing if requested
    if (parallel) {
        bp_param <- MulticoreParam(workers = min(B, 4)) # Limit to 4 cores for Bioconductor
    }

    # Decoupling function for each iteration
    decouple_iteration <- function(b) {
        set.seed(seed + b)
        # Generate independent copy Y from the same distribution
        y_indices <- sample(seq_len(n), size = n, replace = TRUE)
        y <- x[y_indices]

        # Calculate decoupled statistic using C++ implementation
        decoupled_sum <- compute_decoupled_sum(x, y, kernel_function, kernel@symmetric)
        return(decoupled_sum)
    }

    # Run decoupling iterations
    if (parallel) {
        decoupled_samples <- unlist(bplapply(seq_len(B), decouple_iteration, BPPARAM = bp_param))
    } else {
        decoupled_samples <- sapply(seq_len(B), decouple_iteration)
    }

    # Calculate p-value (two-tailed test)
    null_mean <- mean(decoupled_samples)
    null_sd <- sd(decoupled_samples)
    z_score <- (original_stat - null_mean) / null_sd
    p_value <- 2 * (1 - pnorm(abs(z_score)))

    # Create result object
    result <- new("DecoupleResult",
        original_stat = original_stat,
        decoupled_distribution = decoupled_samples,
        kernel_name = kernel_name,
        method = "Friedman-de la Pena Decoupling",
        p_value = p_value
    )

    return(result)
}

#' @useDynLib UStatDecouple, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
