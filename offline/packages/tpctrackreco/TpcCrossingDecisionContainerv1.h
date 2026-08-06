#ifndef TPCTRACKRECO_TPCCROSSINGDECISIONCONTAINERV1_H
#define TPCTRACKRECO_TPCCROSSINGDECISIONCONTAINERV1_H


#include "TpcCrossingDecisionContainer.h"

#include <iostream>
#include <map>

class TpcCrossingDecision;

class TpcCrossingDecisionContainerv1 : public TpcCrossingDecisionContainer
{
 public:
  TpcCrossingDecisionContainerv1();
  ~TpcCrossingDecisionContainerv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override { return static_cast<unsigned int>(m_decisions.size()); }
  void add_decision(TpcCrossingDecision* decision) override;
  TpcCrossingDecision* get_decision(unsigned int assembled_track_id) override;
  const TpcCrossingDecision* get_decision(unsigned int assembled_track_id) const override;

 private:
  std::map<unsigned int, TpcCrossingDecision*> m_decisions;

  ClassDefOverride(TpcCrossingDecisionContainerv1, 1)
};

#endif  // TPCTRACKRECO_TPCCROSSINGDECISIONCONTAINERV1_H
