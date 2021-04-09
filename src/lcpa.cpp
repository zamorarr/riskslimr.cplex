//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <ilcplex/ilocplex.h>
#include "LossComputer.h"
#include "Callback.h"

// [[Rcpp::export]]
Rcpp::List lcpa_cpp(arma::mat x, arma::vec y, arma::vec weights, std::string logfile, int R_max = 3, int time_limit = 60) {
  // create environment
  IloEnv env;

  // initialize the loss computer
  // add intercept to x beforehand
  double c0 = 1E-6;
  LossComputer computer(x, y, weights, c0);

  // inputs
  int d = x.n_cols;
  int R_min = 1;
  R_max = std::min(R_max, d);
  int intercept_min = -10;
  int intercept_max = 10;
  int coef_min = -5;
  int coef_max = 5;
  double L_min = 0;
  double L_max = IloInfinity;
  double V_min = 0;
  double V_max = IloInfinity;

  try {
    // create model
    IloModel model(env);

    // create variables
    IloBoolVarArray alpha(env, d);
    IloIntVarArray lambda(env, d);
    IloIntVar R(env, R_min, R_max);
    IloNumVar L(env, L_min, L_max);
    IloNumVar V(env, V_min, V_max);

    // add constraints
    for (IloInt j = 0; j < d; j++) {

      if (j == 0) {
        lambda[j] = IloIntVar(env, intercept_min, intercept_max);
        model.add(lambda[j] - intercept_min*alpha[j] >= 0);
        model.add(lambda[j] - intercept_max*alpha[j] <= 0);
      } else {
        lambda[j] = IloIntVar(env, coef_min, coef_max);
        model.add(lambda[j] - coef_min*alpha[j] >= 0);
        model.add(lambda[j] - coef_max*alpha[j] <= 0);
      }
    }

    model.add(R == IloSum(alpha));
    model.add(V == L + c0*R);

    // add objective
    model.add(IloMinimize(env, V));

    // define optimizer
    IloCplex cplex(model);

    // add lazy callback
    // We instantiate a FacilityCallback and set the contextMask parameter.
    LossCutCallback cb(lambda, alpha, L, &computer);
    CPXLONG contextMask = 0;
    contextMask |= IloCplex::Callback::Context::Id::Candidate;
    contextMask |= IloCplex::Callback::Context::Id::GlobalProgress;
    cplex.use(&cb, contextMask);

    // cplex parameters
    cplex.setParam(IloCplex::Param::TimeLimit, time_limit); // 60 seconds
    //cplex.setParam(IloCplex::Param::MIP::Display, 0);

    // set output
    std::ofstream outstream(logfile);
    cplex.setOut(outstream);
    //cplex.setOut(Rcpp::Rcout);
    //cplex.setOut(env.getNullStream());
    cplex.setError(Rcpp::Rcerr);
    cplex.setWarning(Rcpp::Rcout);

    // solve model
    if ( !cplex.solve() ) {
      Rcpp::stop("Failed to optimize LP.");
    }

    // write out model for inspection
    //cplex.exportModel("model2.lp");

    // get model values and convert to vectors
    IloNumArray alpha_vals(env);
    cplex.getValues(alpha_vals, alpha);
    std::vector<int> alpha_vec(d);

    IloNumArray lambda_vals(env);
    cplex.getValues(lambda_vals, lambda);
    std::vector<int> lambda_vec(d);

    for (int j = 0; j < d; j++) {
      alpha_vec[j] = alpha_vals[j];
      lambda_vec[j] = lambda_vals[j];
    }

    Rcpp::List result = Rcpp::List::create(
      //Rcpp::Named("status") = cplex.getStatus(),
      Rcpp::Named("objective_value") = cplex.getObjValue(),
      Rcpp::Named("optimality_gap") = cplex.getMIPRelativeGap(),
      Rcpp::Named("alpha") = alpha_vec,
      Rcpp::Named("lambda") = lambda_vec
    );

    env.end();
    return result;
  }
  catch (IloException& e) {
    Rcpp::stop("Concert exception caught");
  }
  catch (...) {
    Rcpp::stop("Unknown exception caught");
  }

  env.end();
  return Rcpp::List::create();
}

