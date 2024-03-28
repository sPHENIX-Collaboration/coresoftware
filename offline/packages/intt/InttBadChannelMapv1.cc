#include "InttBadChannelMapv1.h"

#include <cdbobjects/CDBTTree.h>

#include <string>

InttBadChannelMapv1::~InttBadChannelMapv1()
{
  delete m_bad_channel_set;
}

void InttBadChannelMapv1::identify(
    std::ostream& out) const
{
  out << "InttBadChannelMapv1\n"
      << "\tsize: " << size() << std::endl;
}

std::size_t InttBadChannelMapv1::size() const
{
  if (!m_bad_channel_set)
  {
    return 0;
  }
  return m_bad_channel_set->size();
}

int InttBadChannelMapv1::v_LoadFromCDBTTree(
    CDBTTree& cdbttree)
{
  delete m_bad_channel_set;
  m_bad_channel_set = new Set_t;

  m_bad_channel_set->clear();
  Long64_t N = cdbttree.GetSingleIntValue("size");
  for (Long64_t n = 0; n < N; ++n)
  {
    m_bad_channel_set->insert((struct InttMap::Offline_s){
        .layer =      cdbttree.GetIntValue(n, "layer"),
        .ladder_phi = cdbttree.GetIntValue(n, "ladder_phi"),
        .ladder_z =   cdbttree.GetIntValue(n, "ladder_z"),
        .strip_z =    cdbttree.GetIntValue(n, "strip_z"),
        .strip_phi =  cdbttree.GetIntValue(n, "strip_phi"),
    });
  }

  return 0;
}

bool InttBadChannelMapv1::IsBad(
    InttMap::Offline_s const& ofl) const
{
  if (!m_bad_channel_set)
  {
    return false;
  }
  return m_bad_channel_set->find(ofl) != m_bad_channel_set->end();
}
