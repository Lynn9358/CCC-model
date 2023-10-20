#include <Rcpp.h>
using namespace Rcpp;

// Define the objective function for optimization (SVC_model4)
// This is the same as your previous R function
NumericVector SVC_model4(NumericVector params, List obj, List lr, List recepter, List ligand, double lambda1, double lambda2, double tau, double phi) {
  // ... (Your existing SVC_model4 implementation here)
}
/ BFGS Optimization Function using Rcpp
// [[Rcpp::export]]
NumericVector bfgs_optimize(NumericVector initial_params, List obj, List lr, List recepter, List ligand, double lambda1, double lambda2, double tau, double phi, double tol = 1e-6, int max_iter = 100) {
  int n = initial_params.size();
  NumericVector params = clone(initial_params);

  for (int iter = 0; iter < max_iter; ++iter) {
    NumericVector gradient(n);

    for (int i = 0; i < n; ++i) {
      // Calculate the gradient by finite differences
      double h = 1e-5; // Small step size for finite differences
      NumericVector params_plus_h = clone(params);
      NumericVector params_minus_h = clone(params);
      params_plus_h[i] += h;
      params_minus_h[i] -= h;
      gradient[i] = (SVC_model4(params_plus_h, obj, lr, recepter, ligand, lambda1, lambda2, tau, phi)[0] - SVC_model4(params_minus_h, obj, lr, recepter, ligand, lambda1, lambda2, tau, phi)[0]) / (2 * h);
    }

    if (max(abs(gradient)) < tol) {
      break; // Convergence criterion met
    }

    NumericVector direction(n);

    for (int i = 0; i < n; ++i) {
      direction[i] = -gradient[i];
    }

    // Line search to find the step size
    double step_size = 1.0;
    double c = 0.5; // Backtracking parameter
    double rho = 0.5; // Backtracking parameter

    while (SVC_model4(params + step_size * direction, obj, lr, recepter, ligand, lambda1, lambda2, tau, phi)[0] > SVC_model4(params, obj, lr, recepter, ligand, lambda1, lambda2, tau, phi)[0] + c * step_size * sum(gradient * direction)) {
      step_size *= rho;
    }

    NumericVector s = step_size * direction;

    // Update the parameters
    params += s;
  }

  return params;
}
