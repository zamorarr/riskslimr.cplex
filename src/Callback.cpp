#include "Callback.h"
#include <ilcplex/ilocplex.h>

// see {CPLEX_DIR}/cplex/examples/src/cpp/iloadmipex9.cpp for a heuristic example

// Implementation of the invoke method.
//
// This is the method that we have to implement to fulfill the
// generic callback contract. CPLEX will call this method during
// the solution process at the places that we asked for.
void LossCutCallback::invoke(const IloCplex::Callback::Context &context) {
  if (context.inCandidate()) {
    lazyCut(context);
  } else if (context.inGlobalProgress()) {
    showProgress(context);
  }
}

// show progress callback
void LossCutCallback::showProgress(const IloCplex::Callback::Context &context) {
  time_current = std::time(0);
  if (time_current - time_last_update >= 10) {
    // update every 10 seconds
    time_last_update = time_current;
    IloNum V = context.getIncumbentObjective();

    IloEnv env = context.getEnv();
    //env.setOut(Rcpp::Rcout); // this is causing a C stack overflow error?
    env.out() << "[" << time_current << "] objective value: " << V << "\n";
  }
}


// lazycut callback
void LossCutCallback::lazyCut(const IloCplex::Callback::Context &context) {
  IloEnv env = context.getEnv();

  // discrete coordinate descent on incumbent?
  discreteCoordDescent(context);

  //env.out() << "candidate lambda: [" ;
  arma::vec lambda_val(lambda.getSize());
  for (int j = 0; j < lambda.getSize(); j++) {
    lambda_val[j] = context.getCandidatePoint(lambda[j]);
  }

  // calculate actual loss at this lambda
  // compare that to the candidate objective
  // if actual loss is still greater than candidate loss,
  // add a new cut to improve the candidate
  // us
  //
  // Vmin = candidate? incumbent? value of objective function (uses cutting planes for loss function)
  // Vmax = best value using actual loss function?

  // add lazy constraint
  double loss_point = computer->loss(lambda_val);
  arma::vec loss_slope = computer->loss_grad(lambda_val);

  // check that loss_point is ?
  double L_val = context.getCandidatePoint(L);
  if (abs(loss_point - L_val) <= 1E-8) {
    //env.out() << "close enough, not adding a new constraint" << endl;
    return;
  }

  //  L >= loss_approx + sum_expr(colwise(loss_grad_approx[j]) * lambda[j], j = 1:d) - q
  IloNumExpr sum_expr = IloExpr(env);
  for (int j = 0; j < lambda.getSize(); j++) {
    sum_expr += loss_slope[j]*(lambda[j] - lambda_val[j]);
  }

  context.rejectCandidate(L - loss_point - sum_expr >= 0.0);
  sum_expr.end();
}

void LossCutCallback::discreteCoordDescent(const IloCplex::Callback::Context &context) {
  IloEnv env = context.getEnv();

  // get current candidate values
  IloNum V = context.getCandidateObjective();

  // convert to arma for computations
  IloInt d = lambda.getSize();
  arma::vec lambda_vec(d);
  arma::vec alpha_vec(d);
  arma::vec ub(d);
  arma::vec lb(d);
  for (int j = 0; j < d; j++) {
    lambda_vec[j] = context.getCandidatePoint(lambda[j]);
    alpha_vec[j] = context.getCandidatePoint(alpha[j]);
    ub[j] = lambda[j].getUB();
    lb[j] = lambda[j].getLB();
  }

  arma::vec lambda_vec_best(lambda_vec);
  arma::vec alpha_vec_best(alpha_vec);

  // do one pass of DCD.
  double V_best = V;
  double V_best_new = V_best;

  int iters = 0;
  while(iters < 100) {
    iters++;
    V_best_new = dcdOne(lambda_vec, alpha_vec, lb, ub, lambda_vec_best, alpha_vec_best, V);

    // if no improvement stop looping
    if (V_best_new >= (V_best - 1E-8)) break;

    // update best
    V_best = V_best_new;
    lambda_vec = lambda_vec_best;
    alpha_vec = alpha_vec_best;
    lambda_vec_best.zeros();
    alpha_vec_best.zeros();

  }

  // if V_best is better than V do something
  if (V_best >= (V - 1E-8)) return;

  // convert vec to Iloarray
  IloNumArray lambda_vals(env, d);
  IloNumArray alpha_vals(env, d);
  IloNumArray lambda_vals_best(env, d);
  IloNumArray alpha_vals_best(env, d);
  for (int j = 0; j < d; j++) {
    lambda_vals[j] = context.getCandidatePoint(lambda[j]);
    alpha_vals[j] = context.getCandidatePoint(alpha[j]);
    lambda_vals_best[j] = lambda_vec_best[j];
    alpha_vals_best[j] = alpha_vec_best[j];
  }

  // check now V_best vs V
  //env.out() << "===================" << std::endl;
  //env.out() << "iterations: " << iters << std::endl;
  //env.out() << "candidate objective (cutting planes): " << V << std::endl;
  //env.out() << "candidate lambda: " << lambda_vals << std::endl;
  //env.out() << "candidate alpha: " << alpha_vals << std::endl;
  //env.out() << "dcd objective: " << V_best << std::endl;
  //env.out() << "dcd lambda: " << lambda_vals_best << std::endl;
  //env.out() << "dcd alpha: " << alpha_vals_best << std::endl;

  // post heuristic solution
  IloIntVarArray vars(env);
  vars.add(lambda);
  vars.add(alpha);
  IloNumArray vals(env);
  vals.add(lambda_vals_best);
  vals.add(alpha_vals_best);
  context.postHeuristicSolution(vars, vals, V_best, IloCplex::Callback::Context::SolutionStrategy::Propagate);

  return;
}

double LossCutCallback::dcdOne(
    arma::vec lambda_vec, arma::vec alpha_vec,
    arma::vec lb, arma::vec ub,
    arma::vec &lambda_vec_best, arma::vec &alpha_vec_best,
    double V) {

  // discrete coord descent best
  double V_best = V;
  double V_current = 0;
  int d = lambda_vec.size();
  arma::vec lambda_vec_current(d);
  arma::vec alpha_vec_current(d);

  // for each lambda_j
  for (int j = 0; j < d; j++) {
    // find the list of feasible moves (subject to lb/ub)
    lambda_vec_current = lambda_vec;
    alpha_vec_current = alpha_vec;

    // move down/up along direction j
    if (lambda_vec[j] - 1 >= lb[j]) {
      // set new solution
      lambda_vec_current[j] = lambda_vec[j] - 1;
      alpha_vec_current[j] = lambda_vec_current[j] != 0;

      // compute loss at this solution
      V_current = computer->objective(lambda_vec_current, alpha_vec_current);

      // check if this is better than V_best
      if (V_current < (V_best - 1E-8)) {
        V_best = V_current;
        lambda_vec_best = lambda_vec_current;
        alpha_vec_best = alpha_vec_current;
      }
    }

    // move up
    if (lambda_vec[j] + 1 <= ub[j]) {
      // set new solution
      lambda_vec_current[j] = lambda_vec[j] + 1;
      alpha_vec_current[j] = lambda_vec_current[j] != 0;

      // compute loss at this solution
      V_current = computer->objective(lambda_vec_current, alpha_vec_current);

      // check if this is better than V_best
      if (V_current < (V_best - 1E-8)) {
        V_best = V_current;
        lambda_vec_best = lambda_vec_current;
        alpha_vec_best = alpha_vec_current;
      }
    }
  }

  return V_best;
}
