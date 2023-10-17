#include "MvtxRawRunHeader.h"

void MvtxRawRunHeader::Reset()
{
  m_MvtxL1TrgSet.clear();
}

void MvtxRawRunHeader::identify(std::ostream &os) const
{
  os << "MvtxRawRunHeader" << std::endl;
  os << "Mvtx L1 triggers:  " << m_MvtxL1TrgSet.size()<< std::endl;
}

int MvtxRawRunHeader::isValid() const
{
  return  m_MvtxL1TrgSet.size();
}

void MvtxRawRunHeader::AddL1Trg( const std::set<uint64_t>& mvtxL1TrgSet)
{
  m_MvtxL1TrgSet.insert(mvtxL1TrgSet.cbegin(), mvtxL1TrgSet.cend());
}
