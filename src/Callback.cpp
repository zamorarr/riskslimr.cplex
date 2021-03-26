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

  IloNumArray lambda_vals(env);
  IloNumArray alpha_vals(env);
  context.getCandidatePoint(lambda, lambda_vals);
  context.getCandidatePoint(alpha, alpha_vals);

  // while v < V
  // for each lambda_j
  IloInt d = lambda.getSize();
  for (IloInt j = 0; j < d; j++) {
    // find the list of feasible moves (subject to lb/ub)
    // choose the move that produces the min V
    // end for
  }

  // choose the lambda_j and move that produces the min V
  // update lambda vector to make that move
  // v = V(lambda + delta)

  // output: a new lambda vector with a better V
  // put it postHeuristicSolution

}
