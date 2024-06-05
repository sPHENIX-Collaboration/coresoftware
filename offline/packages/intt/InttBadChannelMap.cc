#include "InttBadChannelMap.h"

#include <cdbobjects/CDBTTree.h>

#include <iostream>
#include <string>

int InttBadChannelMap::LoadFromCdbTTree(CDBTTree& cdbttree)
{
  m_bad_channel_set.clear();
  Long64_t N = cdbttree.GetSingleIntValue("size");
  for (Long64_t n = 0; n < N; ++n)
  {
    m_bad_channel_set.insert((struct InttMap::Offline_s){
        .layer =      cdbttree.GetIntValue(n, "layer"),
        .ladder_phi = cdbttree.GetIntValue(n, "ladder_phi"),
        .ladder_z =   cdbttree.GetIntValue(n, "ladder_z"),
        .strip_z =    cdbttree.GetIntValue(n, "strip_z"),
        .strip_phi =  cdbttree.GetIntValue(n, "strip_phi"),
    });
  }

  return 0;
}

bool InttBadChannelMap::IsBad(InttMap::Offline_s const& ofl) const
{
  return m_bad_channel_set.find(ofl) != m_bad_channel_set.end();
}
