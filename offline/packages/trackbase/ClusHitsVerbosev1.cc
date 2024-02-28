#include "ClusHitsVerbosev1.h"

#include <phool/phool.h>
#include <algorithm>
#include <iostream>

using std::cout;
using std::endl;

namespace
{
  ClusHitsVerbose::Vector dummy_vec;
  ClusHitsVerbose::PairVector dummy_pairvec;
}  // namespace

void ClusHitsVerbosev1::Reset()
{
  m_data.clear();
  m_stage_phi.clear();
  m_stage_z.clear();
  m_stage_phiCut.clear();
  m_stage_zCut.clear();
}

bool ClusHitsVerbosev1::hasClusKey(TrkrDefs::cluskey key) const
{
  return m_data.find(key) != m_data.end();
}

ClusHitsVerbose::Vector& ClusHitsVerbosev1::vecBins(TrkrDefs::cluskey key, int which)
{
  if (!hasClusKey(key))
  {
    return dummy_vec;
  }
  else
  {
    return m_data[key][which];
  }
}

ClusHitsVerbose::Vector& ClusHitsVerbosev1::phiBins(TrkrDefs::cluskey key)
{
  return vecBins(key, 0);
}

ClusHitsVerbose::Vector& ClusHitsVerbosev1::zBins(TrkrDefs::cluskey key)
{
  return vecBins(key, 1);
}

ClusHitsVerbose::Vector& ClusHitsVerbosev1::phiCutBins(TrkrDefs::cluskey key)
{
  return vecBins(key, 2);
}

ClusHitsVerbose::Vector& ClusHitsVerbosev1::zCutBins(TrkrDefs::cluskey key)
{
  return vecBins(key, 3);
}

ClusHitsVerbose::PairVector ClusHitsVerbosev1::pvecIE(TrkrDefs::cluskey key, int which)
{
  if (!hasClusKey(key))
  {
    return dummy_pairvec;
  }

  // unzip the vector<pair<int,int>> into two vectors: vector<int> vector<int>
  std::vector<int> vecI{};
  std::vector<int> vecE{};
  for (auto& entry : m_data[key][which])
  {
    vecI.push_back(entry.first);
    vecE.push_back(entry.second);
  }
  return {vecI, vecE};
}

ClusHitsVerbose::PairVector ClusHitsVerbosev1::phiBins_pvecIE(TrkrDefs::cluskey key)
{
  return pvecIE(key, 0);
}

ClusHitsVerbose::PairVector ClusHitsVerbosev1::zBins_pvecIE(TrkrDefs::cluskey key)
{
  return pvecIE(key, 1);
}

ClusHitsVerbose::PairVector ClusHitsVerbosev1::phiCutBins_pvecIE(TrkrDefs::cluskey key)
{
  return pvecIE(key, 2);
}

ClusHitsVerbose::PairVector ClusHitsVerbosev1::zCutBins_pvecIE(TrkrDefs::cluskey key)
{
  return pvecIE(key, 3);
}

void ClusHitsVerbosev1::push_hits(TrkrDefs::cluskey key)
{
  m_data[key] = {m_stage_phi, m_stage_z, m_stage_phiCut, m_stage_zCut};
  m_stage_phi.clear();
  m_stage_z.clear();
  m_stage_phiCut.clear();
  m_stage_zCut.clear();
}
