#include "ClusHitsVerbose.h"

namespace
{
  ClusHitsVerbose::Vector dummy_vec;
  ClusHitsVerbose::Map dummy_map;
  ClusHitsVerbose::PairVector dummy_pairvec;
}  // namespace

ClusHitsVerbose::Vector& ClusHitsVerbose::phiBins(TrkrDefs::cluskey /*this key*/)
{
  return dummy_vec;
}

ClusHitsVerbose::Vector& ClusHitsVerbose::zBins(TrkrDefs::cluskey /*this key*/)
{
  return dummy_vec;
}

ClusHitsVerbose::Vector& ClusHitsVerbose::phiCutBins(TrkrDefs::cluskey /*this key*/)
{
  return dummy_vec;
}

ClusHitsVerbose::Vector& ClusHitsVerbose::zCutBins(TrkrDefs::cluskey /*this key*/)
{
  return dummy_vec;
}

ClusHitsVerbose::Map& ClusHitsVerbose::getMap()
{
  return dummy_map;
}

ClusHitsVerbose::PairVector ClusHitsVerbose::phiBins_pvecIE(TrkrDefs::cluskey /*this key*/)
{
  return dummy_pairvec;
}
ClusHitsVerbose::PairVector ClusHitsVerbose::phiCutBins_pvecIE(TrkrDefs::cluskey /*this key*/)
{
  return dummy_pairvec;
}
ClusHitsVerbose::PairVector ClusHitsVerbose::zBins_pvecIE(TrkrDefs::cluskey /*this key*/)
{
  return dummy_pairvec;
}
ClusHitsVerbose::PairVector ClusHitsVerbose::zCutBins_pvecIE(TrkrDefs::cluskey /*this key*/)
{
  return dummy_pairvec;
}
