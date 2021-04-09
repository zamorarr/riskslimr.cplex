#include "LossComputer.h"

LossComputer::LossComputer(const arma::mat &_x, const arma::vec &_y, const arma::vec &_w, double _c0):
  x(_x), y(_y), w(_w), c0(_c0) {

  z = x.each_col() % y;
  n = x.n_rows;
  d = x.n_cols;
  w_sum = arma::sum(w);
}

double LossComputer::loss(arma::vec lambda) {
  //return arma::mean(arma::log(1 + arma::exp(-z * lambda)));
  arma::vec loss_i = arma::log(1 + arma::exp(-z * lambda));
  return arma::sum(w % loss_i)/w_sum;
}

arma::vec LossComputer::loss_grad(arma::vec lambda) {
  arma::vec b = 1 + arma::exp(z * lambda); // {n x 1} vector
  //arma::mat a = z.each_col() / b; // {n x d} matrix
  //arma::mat lg = arma::mean(-a, 0); // col means {1 x d} matrix
  arma::mat a = z.each_col() % (w / b);
  arma::mat lg = -arma::sum(a, 0)/w_sum; // col sums
  return arma::vectorise(lg);
}

double LossComputer::objective(arma::vec lambda, arma::vec alpha) {
  int R = arma::sum(alpha);
  return loss(lambda) + c0*R;
}

// [[Rcpp::export]]
double test_loss_cpp(arma::mat x, arma::vec y, arma::vec w, double c0, arma::vec lambda) {
  LossComputer computer(x, y, w, c0);
  return computer.loss(lambda);
}

// [[Rcpp::export]]
arma::vec test_loss_grad_cpp(arma::mat x, arma::vec y, arma::vec w, double c0, arma::vec lambda) {
  LossComputer computer(x, y, w, c0);
  return computer.loss_grad(lambda);
}

