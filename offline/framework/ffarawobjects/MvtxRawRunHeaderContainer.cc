#include "MvtxRawRunHeaderContainer.h"

void MvtxRawRunHeaderContainer::Reset()
{
  m_MvtxL1TrgSet.clear();
}

void MvtxRawRunHeaderContainer::identify(std::ostream &os) const
{
  os << "MvtxRawRunHeaderContainer" << std::endl;
  os << "Mvtx L1 triggers:  " << m_MvtxL1TrgSet.size()<< std::endl;
}

int MvtxRawRunHeaderContainer::isValid() const
{
  return  m_MvtxL1TrgSet.size();
}

void MvtxRawRunHeaderContainer::AddL1Trg( const std::set<uint64_t>& mvtxL1TrgSet)
{
  m_MvtxL1TrgSet.insert(mvtxL1TrgSet.cbegin(), mvtxL1TrgSet.cend());
}
