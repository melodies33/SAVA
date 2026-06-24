#include <Rcpp.h>
using namespace Rcpp;

// Function: gamldpp
//
// Computes the gamma sequence used by the LORD++ procedure (Ramdas et al., 2017).
//
// Inputs:
//   x      : index of gamma.
//
// Returns:
//   The gamma value of index x.
//
// See Appendix B.2 for details of the gamma sequence definition.

double gamldpp(double x) {
  double rr = 0.0722*log(std::max(x,2.0))/x/exp(sqrt(log(x)));
  return rr;
}

// Function: ldpp
//
// Computes the test levels in the LORD++  procedure (Ramdas et al., 2017) given the index vector of gamma.
//
// Inputs:
//   q      : target FDR level.
//   W      : initial alpha wealth.
//   i      : current task index.
//   lags   : the index vector of gamma.
//
// Returns:
//   The test level for the current task used in LORD++ given the index vector of gamma.
//
// See Appendix B.2 for details of the calculation.

double ldpp(double q, double W, double i, NumericVector lags){
  double re = W*gamldpp(i);
  int ll = lags.size();
  if (ll > 0){
    re += (q - W)*gamldpp(lags[0]);
    if (ll > 1){
      for (int l = 1; l < ll; l++){
        re += q * gamldpp(lags[l]);
      }
    }
  }
  return re;
}

// Function: threldpp
//
// Computes the LORD++ test levels used in the simulation comparisons with the proposed SAVA procedure.
//
// Inputs:
//   q      : target FDR level.
//   W      : initial alpha wealth.
//   index  : current task index.
//   rej    : vector of previous rejection indices.
//
// Returns:
//   The LORD++ test level for the current task.
//
// See Appendix B.2 for details.

//[[Rcpp::export]]
double threldpp(double q, double W, double index, NumericVector rej){
  NumericVector cc = rej[rej < index];
  cc = index - cc;
  return ldpp(q, W, index, cc);
}

// Function: gamsaff
//
// Computes the gamma sequence used by the SAFFRON (Ramdas et al., 2018) and ADDIS (Tian and Ramdas, 2019) procedures.
//
// Inputs:
//   i      : index of gamma.
//
// Returns:
//   The gamma value of index i.
//
// See Appendix B.2 for details.

double gamsaff(double i){
  return 0.4374901658/pow(i, 1.6);
}


// Function: sffr
//
// Computes the test levels for the SAFFRON procedure (Ramdas et al. 2018) given the vector of rejection indices. 
//
// Inputs:
//   q      : target FDR level.
//   lamb   : parameter lambda used in SAFFRON.
//   w0     : initial alpha wealth.
//   lags   : the index vector of gamma.
//   index  : current task index.
//
// Returns:
//   The SAFFRON test level for the current task given the indices of gamma.
//
// See Appendix B.2 for details of the calculation.

double sffr(double q, double lamb, double w0, NumericVector lags, double index){
  double l = lags.length();
  double re = w0*gamsaff(lags[0]);
  if (l > 1){
    re += (q - w0)*gamsaff(lags[1]);
    if (l > 2){
      for (int ll = 2; ll < l; ll++){
        re += q*gamsaff(lags[ll]);
      }
    }
  }
  re = std::min(lamb, re*(1 - lamb));
  return re;
}


// Function: csaff
//
// Computes the sum of tau_j and C_{j+} required by the SAFFRON (Ramdas et al., 2018) and ADDIS (Tian and Ramdas, 2019) procedures.
//
// Inputs:
//   index  : current task index.
//   cvec   : vector of candidate indicators C_j.
//   tau    : vector of rejection indices.
//   
//
// Returns:
//   The vector of the sum of tau_j and C_{j+} used in SAFFRON and ADDIS given the current index.
//
// See Appendix B.2 for details.

NumericVector csaff(double index, NumericVector cvec, NumericVector tau){
  int taulen = tau.length();
  NumericVector ree(1 + taulen);
  for (int i = 0; i < index - 1; i++){
    ree[0] += cvec[i];
  }
  if (taulen > 0){
    for (int j = 0; j < taulen; j++){
      if (tau[j] + 1 < index){
        for (int k = tau[j]; k < index - 1; k++){
          ree[j+1] += cvec[k];
        }
      }else{
        ree[j+1] = 0;
      }
    }
  }
  return ree;
}

// Function: thresaff
//
// Computes the SAFFRON test levels used in the simulation comparisons with the proposed SAVA procedure.
//
// Inputs:
//   q      : target FDR level.
//   w0     : initial alpha wealth.
//   index  : current task index.
//   rejvec : vector of the indexes of rejections. 
//   cvec   : vector of C_j in SAFFRON.
//   lamb   : parameter lambda used in SAFFRON.
//
// Returns:
//   The SAFFRON test level for the current task.
//
// See Appendix B.2 for details.

//[[Rcpp::export]]
double thresaff(double q, double w0, double index, NumericVector rejvec, NumericVector cvec, double lamb){
  NumericVector cc = rejvec[rejvec < index];
  NumericVector cff = csaff(index, cvec, cc);
  NumericVector rejlag = index - cff;
  double thre;
  int ll = cc.length();
  if (ll > 0){
    for (int k = 1; k < (ll + 1); k++){
      rejlag[k] = rejlag[k] - cc[k - 1];
    }
  }
  thre= sffr(q, lamb, w0, rejlag, index);
  return thre;
}

// Function: threaddis
//
// Computes the ADDIS test levels used in the simulation comparisons with the proposed SAVA procedure.
//
// Inputs:
//   q      : target FDR level.
//   w0     : initial alpha wealth.
//   index  : current task index.
//   rejvec : vector of the indexes of rejections. 
//   cvec   : vector of C_j in ADDIS.
//   sindex : count S_j in ADDIS.
//   kstar  : vector of kappa_i^* quantities used by ADDIS.
//   lamb   : parameter lambda used in ADDIS.
//   tau    : parameter tau used in ADDIS.
//
// Returns:
//   The ADDIS test level for the current task.
//
// See Appendix B.2 for details.

//[[Rcpp::export]]
double threaddis(double q, double w0, double index, NumericVector rejvec, NumericVector cvec, double sindex, NumericVector kstar, double lamb, double tau){
  NumericVector cc = rejvec[rejvec < index];
  NumericVector cff = csaff(index, cvec, cc);
  NumericVector rejlag = sindex - cff;
  int ll = cc.length();
  if (ll > 0){
    for (int k = 1; k < (ll + 1); k++){
      rejlag[k] = rejlag[k] - kstar[k - 1];
    }
  }
  double l = rejlag.length();
  double re = w0*gamsaff(rejlag[0] + 1);
  if (l > 1){
    re += (q - w0)*gamsaff(rejlag[1] + 1);
    if (l > 2){
      for (int ll = 2; ll < l; ll++){
        re += q*gamsaff(rejlag[ll] + 1);
      }
    }
  }
  re = std::min(lamb, re*(tau - lamb));
  return re;
}

