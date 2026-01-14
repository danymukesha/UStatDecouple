#' @title DecoupleResult Class
#' @description S4 class storing results from decoupling analysis
#' @slot original_stat Numeric value of the original U-statistic
#' @slot decoupled_distribution Numeric vector of decoupled statistics
#' @slot kernel_name Character string describing the kernel used
#' @slot method Character string describing the decoupling method
#' @slot p_value Numeric p-value from hypothesis testing
#' @export
setClass("DecoupleResult",
    slots = list(
        original_stat = "numeric",
        decoupled_distribution = "numeric",
        kernel_name = "character",
        method = "character",
        p_value = "numeric"
    )
)

#' @title UStatKernel Class
#' @description S4 class for U-statistic kernel functions
#' @slot kernel_function Function implementing the kernel
#' @slot kernel_name Character string name of the kernel
#' @slot symmetric Logical indicating if kernel is symmetric
#' @export
setClass("UStatKernel",
    slots = list(
        kernel_function = "function",
        kernel_name = "character",
        symmetric = "logical"
    )
)

#' @title Show method for DecoupleResult
#' @param object DecoupleResult object to display
#' @export
setMethod("show", "DecoupleResult", function(object) {
    cat("DecoupleResult object:\n")
    cat(sprintf("  Original U-statistic: %.4f\n", object@original_stat))
    cat(sprintf("  Decoupled mean: %.4f\n", mean(object@decoupled_distribution)))
    cat(sprintf("  Decoupled SD: %.4f\n", sd(object@decoupled_distribution)))
    cat(sprintf("  Kernel: %s\n", object@kernel_name))
    cat(sprintf("  Method: %s\n", object@method))
    cat(sprintf("  P-value: %.4f\n", object@p_value))

    # Calculate z-score and significance
    z_score <- (object@original_stat - mean(object@decoupled_distribution)) /
        sd(object@decoupled_distribution)
    cat(sprintf("  Z-score: %.4f\n", z_score))

    if (!is.na(object@p_value)) {
        significance <- ifelse(object@p_value < 0.001, "***",
            ifelse(object@p_value < 0.01, "**",
                ifelse(object@p_value < 0.05, "*", " ")
            )
        )
        cat(sprintf("  Significance: %s (p = %.4f)\n", significance, object@p_value))
    }
})

#' @title Create U-statistic kernel
#' @param kernel_function Function taking two arguments
#' @param kernel_name Character name for the kernel
#' @param symmetric Logical indicating if kernel is symmetric
#' @return UStatKernel object
#' @export
create_kernel <- function(kernel_function, kernel_name, symmetric = TRUE) {
    new("UStatKernel",
        kernel_function = kernel_function,
        kernel_name = kernel_name,
        symmetric = symmetric
    )
}
