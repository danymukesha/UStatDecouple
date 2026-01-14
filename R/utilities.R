#' @title Load Example DNA Sequences
#' @description Load example DNA sequences from package data
#' @return List of DNA sequences
#' @export
load_example_sequences <- function() {
    # This would typically load real data from inst/extdata/
    # For demonstration, we'll create a small example
    bases <- c("A", "C", "G", "T")
    set.seed(42)

    sequences <- list(
        seq1 = c("A", "T", "G", "C", "A", "T", "G", "C"),
        seq2 = c("A", "T", "G", "T", "A", "T", "G", "C"),
        seq3 = c("A", "C", "G", "C", "A", "T", "G", "T"),
        seq4 = c("A", "T", "C", "C", "A", "T", "G", "C"),
        seq5 = c("A", "T", "G", "C", "G", "T", "G", "C")
    )

    return(sequences)
}

#' @title Validate Input Data
#' @description Validate input data for decoupling analysis
#' @param x Input data to validate
#' @param kernel Kernel to validate against
#' @export
validate_input_data <- function(x, kernel) {
    if (!is.list(x)) {
        stop("Input data must be a list")
    }

    if (length(x) < 2) {
        stop("Sample size must be at least 2")
    }

    # Check that all elements have the same structure
    element_types <- sapply(x, class)
    if (length(unique(element_types)) > 1) {
        warning("Input elements have different types/classes")
    }

    # Validate kernel function
    if (!is.function(kernel@kernel_function)) {
        stop("Kernel function must be a function")
    }

    tryCatch(
        {
            # Test kernel on first two elements
            test_result <- kernel@kernel_function(x[[1]], x[[2]])
            if (!is.numeric(test_result) || length(test_result) != 1) {
                stop("Kernel function must return a single numeric value")
            }
        },
        error = function(e) {
            stop(sprintf("Kernel validation failed: %s", e$message))
        }
    )

    return(TRUE)
}
