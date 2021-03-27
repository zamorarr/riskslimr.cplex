#include "LossComputer.h"

LossComputer::LossComputer(const arma::mat &_x, const arma::vec &_y, double _c0):
  x(_x), y(_y), c0(_c0) {
  z = x.each_col() % y;
}

double LossComputer::loss(arma::vec lambda) {
  return arma::mean(arma::log(1 + arma::exp(-z * lambda)));
}

arma::vec LossComputer::loss_grad(arma::vec lambda) {
  arma::vec b = 1 + arma::exp(z * lambda); // {n x 1} vector
  arma::mat a = z.each_col() / b;
  arma::mat lg = arma::mean(-a, 0); // col means
  return arma::vectorise(lg);
}

double LossComputer::objective(arma::vec lambda, arma::vec alpha) {
  int R = arma::sum(alpha);
  return loss(lambda) + c0*R;
}
