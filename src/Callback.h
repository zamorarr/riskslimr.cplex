#ifndef CALLBACK_H
#define CALLBACK_H

#include <ilcplex/ilocplex.h>
#include "LossComputer.h"

// my callback implements the Callback::Function interface (I think?)
class LossCutCallback: public IloCplex::Callback::Function {
private:
  // Empty constructor is forbidden
  LossCutCallback() {};

  // Copy constructor is forbidden
  LossCutCallback(const LossCutCallback &tocopy) {};

  // current solution
  IloIntVarArray lambda;
  IloIntVarArray alpha;
  IloNumVar L;

  // loss computer
  LossComputer* computer;

  // keep track of display information
  std::time_t time_last_update = std::time(0);
  std::time_t time_current = std::time(0);

  // callbacks
  void lazyCut(const IloCplex::Callback::Context &context);
  void discreteCoordDescent(const IloCplex::Callback::Context &context);
  void showProgress(const IloCplex::Callback::Context &context);

public:
  // Constructor with data
  LossCutCallback(const IloIntVarArray &_lambda, IloIntVarArray &_alpha, IloNumVar &_L, LossComputer* _computer):
  lambda(_lambda), alpha(_alpha), L(_L), computer(_computer) {}

  // This is the function that we have to implement and that CPLEX will call
  // during the solution process at the places that we asked for.
  virtual void invoke(const IloCplex::Callback::Context &context) ILO_OVERRIDE;

  // Destructor
  ~LossCutCallback() {};
};

#endif /* CALLBACK_H */
