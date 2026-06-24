#include <Rcpp.h>
using namespace Rcpp;


// Function: lambda
//
// Computes the parameter lambda used in the e-processes implemented in Section 4.1, following the betting-based e-process framework of Waudby-Smith & Ramdas (2023).
//
// Inputs:
//   alpha   : the test level.
//   t       : index of lambda.
//
// Returns:
//   The value of lambda for index t.
//
// See Section 4 for details.

double lambda(double alpha, double t){
  double re = 8*log(2.0/alpha)/(t * log(t + 1));
  re = std::min(1.0, sqrt(re));
  return re;
}


// Function: oracsprt
//
// Computes the e-processes implemented in Section B.3 of Appendix.
//
// Inputs:
//   realmu  : alternative mean.
//   muzero  : null mean.
//   meansamp: sample mean.
//   n       : sample size.
//
// Returns:
//   The value of the oracle e-process (likelihood ratio) comparing muzero against realmu.
//
// See Section B.3.1 for details.

// [[Rcpp::export]]
double oracsprt(double realmu, double muzero, double meansamp, double n){
  double re = exp(n*(realmu - muzero)*(2*meansamp - muzero - realmu)/2);
  return re;
}

// Function: falphaplus
//
// Computes the lower confidence-bound statistic corresponding to a specified confidence level under the betting-based e-process framework of Waudby-Smith and Ramdas (2023).
//
// By Proposition 1 of Waudby-Smith and Ramdas (2023), falphaplus(alpha, ...) > 0 is equivalent to the corresponding e-process for arm A crossing the significance threshold alpha.
//
// This quantity is used for decision making in the simulation studies of Section 4 and Appendix B.4.
//
// Inputs:
//   alpha   : the test level.
//   sampvec : vector of observed samples.
//   bound   : the support width.
//   t       : current sample size.
//   q       : target FSR level.
//
// Returns:
//   The value of the lower confidence-bound statistic associated with the betting-based e-process.
//
// See Section 4 for details.


// [[Rcpp::export]]
double falphaplus(double alpha, NumericVector sampvec, double bound, double t, double q){
  NumericVector lamvec(t);
  double re;
  for (int i = 0; i < t; i++){
    lamvec[i] = lambda(q, (double)(i + 1));
  }
  re = sum(lamvec * sampvec) - (log(2.0/alpha) - sum(pow(lamvec, 2))/8.0)*bound;
  return re;
}

// Function: falphaminus
//
// Computes the upper confidence-bound statistic corresponding to a specified confidence level under the betting-based e-process framework of Waudby-Smith and Ramdas (2023).
//
// By Proposition 1 of Waudby-Smith and Ramdas (2023), falphaminus(alpha, ...) < 0 is equivalent to the corresponding e-process for arm B crossing the significance threshold alpha.
//
// This quantity is used for decision making in the simulation studies of Section 4 and Appendix B.4.
//
// Inputs:
//   alpha   : the test level.
//   sampvec : vector of observed samples.
//   bound   : the support width.
//   t       : current sample size.
//   q       : target FSR level.
//
// Returns:
//   The value of the lower confidence-bound statistic associated with the betting-based e-process.
//
// See Section 4 for details.

// [[Rcpp::export]]
double falphaminus(double alpha, NumericVector sampvec, double bound, double t, double q){
  NumericVector lamvec(t);
  double re;
  for (int i = 0; i < t; i++){
    lamvec[i] = lambda(q, (double)(i + 1));
  }
  re = sum(lamvec * sampvec) + (log(2.0/alpha) + sum(pow(lamvec, 2))/8.0)*bound;
  return re;
}
