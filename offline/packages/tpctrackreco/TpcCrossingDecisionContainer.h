#ifndef TPCTRACKRECO_TPCCROSSINGDECISIONCONTAINER_H
#define TPCTRACKRECO_TPCCROSSINGDECISIONCONTAINER_H


#include <phool/PHObject.h>

#include <iostream>

class TpcCrossingDecision;

class TpcCrossingDecisionContainer : public PHObject
{
 public:
  TpcCrossingDecisionContainer() = default;
  ~TpcCrossingDecisionContainer() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "TpcCrossingDecisionContainer base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  virtual unsigned int size() const { return 0; }
  virtual void add_decision(TpcCrossingDecision*) {}
  virtual TpcCrossingDecision* get_decision(unsigned int) { return nullptr; }
  virtual const TpcCrossingDecision* get_decision(unsigned int) const { return nullptr; }

 private:
  ClassDefOverride(TpcCrossingDecisionContainer, 0)
};

#endif  // TPCTRACKRECO_TPCCROSSINGDECISIONCONTAINER_H
