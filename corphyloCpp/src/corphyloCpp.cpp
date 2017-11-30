#include <RcppArmadillo.h>
#include <numeric>
#include <cmath>

using namespace Rcpp;


typedef arma::uword uint;
typedef arma::sword sint;

#define COND_MIN 0.0000000001
#define MAX_RETURN 10000000000


// Needed for multiple functions below (it's faster than using std::pow(a, b))
// a^b
double fast_pow(const double& a, const double& b) {
    return std::exp(b * std::log(a));
}
arma::vec fast_pow(const double& a, const arma::vec& b) {
    uint n = b.n_elem;
    arma::vec x(n);
    for (uint i = 0; i < n; i++) x(i) = std::exp(b(i) * std::log(a));
    return x;
}
arma::rowvec fast_pow(const double& a, const arma::rowvec& b) {
    uint n = b.n_elem;
    arma::rowvec x(n);
    for (uint i = 0; i < n; i++) x(i) = std::exp(b(i) * std::log(a));
    return x;
}
arma::vec fast_pow(const arma::vec& a, const double& b) {
    uint n = a.n_elem;
    arma::vec x(n);
    for (uint i = 0; i < n; i++) x(i) = std::exp(b * std::log(a(i)));
    return x;
}
arma::mat fast_pow(const double& a, const arma::mat& b) {
    uint nr = b.n_rows, nc = b.n_cols;
    arma::mat x(nr, nc);
    for (uint i = 0; i < nr; i++) {
        for (uint j = 0; j < nc; j++) {
            x(i, j) = std::exp(b(i, j) * std::log(a));
        }
    }
    return x;
}

// Transpose functions that return what I want them to...
inline arma::mat tp(const arma::mat& M){
    return M.t();
}
inline arma::rowvec tp(const arma::vec& V){
    return arma::conv_to<arma::rowvec>::from(V.t());
}
inline arma::vec tp(const arma::rowvec& V){
    return arma::conv_to<arma::vec>::from(V.t());
}



//[[Rcpp::export]]
double corphylo_LL(const arma::vec& par, const arma::mat& XX, const arma::mat& UU, 
                   const arma::mat& MM, const arma::mat& tau, const arma::mat& Vphy, 
                   bool REML = true, bool constrain_d = false, bool verbose = false) {
    
    uint n = Vphy.n_rows;
    uint p = XX.n_rows / n;
    arma::vec L_elements = par(arma::span(0, (p + p * (p - 1)/2) - 1));
    arma::mat L(p, p, arma::fill::zeros);
    for (uint i = 0, j = 0, k = p - 1; i < p; i++) {
        L(arma::span(i, p-1), i) = L_elements(arma::span(j, k));
        j = k + 1;
        k += (p - i - 1);
    }
    arma::mat R = L.t() * L;
    arma::vec d;
    if (constrain_d) {
        arma::vec logit_d = par(arma::span((p + p * (p - 1) / 2), par.n_elem - 1));
        if (arma::max(arma::abs(logit_d)) > 10)  return MAX_RETURN;
        d = 1/(1 + arma::exp(-logit_d));
    } else {
        d = par(arma::span((p + p * (p - 1) / 2), par.n_elem - 1));
        if (max(d) > 10) return MAX_RETURN;
    }
    
    // OU transform
    arma::mat C(p * n, p * n, arma::fill::zeros);
    for (uint i = 0; i < p; i++) {
        arma::mat Cd;
        for (uint j = 0; j < p; j++) {
            Cd = fast_pow(d(i), tau) % fast_pow(d(j), tp(tau)) % 
                (1 - fast_pow(d(i) * d(j), Vphy));
            Cd /= (1 - d(i) * d(j));
            C(arma::span(n * i, (i + 1) * n - 1), arma::span(n * j, (j + 1) * n - 1)) =
                R(i, j) * Cd;
        }
    }
    
    arma::mat V = C;
    V += arma::diagmat(arma::vectorise(MM));
    double rcond_dbl = arma::rcond(V);
    if (!arma::is_finite(rcond_dbl) || rcond_dbl < COND_MIN) return MAX_RETURN;
    
    arma::mat iV = arma::inv(V);
    arma::mat denom = tp(UU) * iV * UU;
    rcond_dbl = arma::rcond(denom);
    if (!arma::is_finite(rcond_dbl) || rcond_dbl < COND_MIN) return MAX_RETURN;
    
    arma::mat num = tp(UU) * iV * XX;
    arma::mat B = arma::solve(denom, num);
    arma::mat H = XX - UU * B;
    
    double logdetV, det_sign;
    arma::log_det(logdetV, det_sign, iV);
    if (!arma::is_finite(logdetV)) return MAX_RETURN;
    logdetV *= -1;
    
    double LL;
    if (REML) {
        arma::mat to_det = tp(UU) * iV * UU;
        double det_val;
        arma::log_det(det_val, det_sign, to_det);
        double lhs = arma::as_scalar(tp(H) * iV * H);
        LL = 0.5 * (logdetV + det_val + lhs);
    } else {
        LL = 0.5 * arma::as_scalar(logdetV + tp(H) * iV * H);
    }
    
    return LL;
}