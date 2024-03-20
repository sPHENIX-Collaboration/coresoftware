#include "InttMaskedChannelSetv1.h"

#include <cdbobjects/CDBTTree.h>

#include <string>

InttMaskedChannelSetv1::~InttMaskedChannelSetv1()
{
  delete m_HotChannelSet;
}

void InttMaskedChannelSetv1::identify(
    std::ostream& out) const
{
  out << "InttMaskedChannelSetv1"
      << "\n"
      << "\tsize: " << size() << std::endl;
}

std::size_t InttMaskedChannelSetv1::size() const
{
  if (!m_HotChannelSet)
  {
    return 0;
  }
  return m_HotChannelSet->size();
}

int InttMaskedChannelSetv1::v_LoadFromCDBTTree(
    CDBTTree& cdbttree)
{
  delete m_HotChannelSet;
  m_HotChannelSet = new Set_t;

  m_HotChannelSet->clear();
  Long64_t N = cdbttree.GetSingleIntValue("size");
  for (Long64_t n = 0; n < N; ++n)
  {
    m_HotChannelSet->insert((struct InttMap::Offline_s){
        .layer = cdbttree.GetIntValue(n, "layer"),
        .ladder_phi = cdbttree.GetIntValue(n, "ladder_phi"),
        .ladder_z = cdbttree.GetIntValue(n, "ladder_z"),
        .strip_z = cdbttree.GetIntValue(n, "strip_z"),
        .strip_phi = cdbttree.GetIntValue(n, "strip_phi"),
    });
  }

  return 0;
}

bool InttMaskedChannelSetv1::IsDeadChannel(
    InttMap::Offline_s const& ofl) const
{
  if (!m_HotChannelSet)
  {
    return false;
  }
  return m_HotChannelSet->find(ofl) != m_HotChannelSet->end();
}
