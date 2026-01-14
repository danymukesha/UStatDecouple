#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double compute_decoupled_sum(List x, List y, Function kernel_function, bool symmetric) {
    int n = x.size();
    double total = 0.0;
    int count = 0;

    if (symmetric) {
        // For symmetric kernels, we only need to compute i < j
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                // Call R kernel function
                double value = as<double>(kernel_function(x[i], y[j]));
                total += value;
                count++;
            }
        }
        // Since we only computed i < j, we need to multiply by 2 for the full symmetric sum
        // But divide by the number of unique pairs
        return total / count;
    } else {
        // For asymmetric kernels, compute all pairs where i != j
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    double value = as<double>(kernel_function(x[i], y[j]));
                    total += value;
                    count++;
                }
            }
        }
        return total / count;
    }
}

// [[Rcpp::export]]
NumericVector compute_multiple_decoupled_sums(List x, List y_list, Function kernel_function, bool symmetric) {
    int B = y_list.size();
    NumericVector results(B);

    for (int b = 0; b < B; b++) {
        List y = y_list[b];
        results[b] = compute_decoupled_sum(x, y, kernel_function, symmetric);
    }

    return results;
}

// You can enable OpenMP by setting PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS) in src/Makevars
// And PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS) in src/Makevars
