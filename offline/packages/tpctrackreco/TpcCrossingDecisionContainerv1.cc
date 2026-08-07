#include "TpcCrossingDecisionContainerv1.h"

#include "TpcCrossingDecision.h"

ClassImp(TpcCrossingDecisionContainerv1)

TpcCrossingDecisionContainerv1::TpcCrossingDecisionContainerv1()
{
  Reset();
}

TpcCrossingDecisionContainerv1::~TpcCrossingDecisionContainerv1()
{
  Reset();
}

void TpcCrossingDecisionContainerv1::identify(std::ostream& os) const
{
  os << "TpcCrossingDecisionContainerv1 with "
     << m_decisions.size()
     << " decisions"
     << std::endl;
}

void TpcCrossingDecisionContainerv1::Reset()
{
  for (auto& decision : m_decisions) { delete decision.second;
}
  m_decisions.clear();
}

int TpcCrossingDecisionContainerv1::isValid() const
{
  return m_decisions.empty() ? 0 : 1;
}

PHObject* TpcCrossingDecisionContainerv1::CloneMe() const
{
  TpcCrossingDecisionContainerv1* copy = new TpcCrossingDecisionContainerv1();
  for (const auto& decision : m_decisions)
  {
    if (decision.second) { copy->m_decisions[decision.first] = static_cast<TpcCrossingDecision*>(decision.second->CloneMe());
}
  }
  return copy;
}

void TpcCrossingDecisionContainerv1::add_decision(TpcCrossingDecision* decision)
{
  if (!decision) { return;
}

  const unsigned int assembled_track_id = decision->get_assembled_track_id();
  auto iter = m_decisions.find(assembled_track_id);
  if (iter != m_decisions.end())
  {
    delete iter->second;
    iter->second = decision;
    return;
  }
  m_decisions[assembled_track_id] = decision;
}

TpcCrossingDecision* TpcCrossingDecisionContainerv1::get_decision(unsigned int assembled_track_id)
{
  auto iter = m_decisions.find(assembled_track_id);
  return iter == m_decisions.end() ? nullptr : iter->second;
}

const TpcCrossingDecision* TpcCrossingDecisionContainerv1::get_decision(unsigned int assembled_track_id) const
{
  auto iter = m_decisions.find(assembled_track_id);
  return iter == m_decisions.end() ? nullptr : iter->second;
}
