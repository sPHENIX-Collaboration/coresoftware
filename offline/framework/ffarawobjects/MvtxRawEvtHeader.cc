#include "MvtxRawEvtHeader.h"

void MvtxRawEvtHeader::Reset()
{
  m_MvtxFeeIdSet.clear();
  m_MvtxL1TrgSet.clear();
}

void MvtxRawEvtHeader::identify(std::ostream &os) const
{
  os << "MvtxRawEvtHeader" << std::endl;
  os << "Feeids in the event:  " << m_MvtxFeeIdSet.size()<< std::endl;
  os << "Mvtx L1 triggers in the event:  " << m_MvtxL1TrgSet.size()<< std::endl;
}

int MvtxRawEvtHeader::isValid() const
{
  return  m_MvtxFeeIdSet.size();
}

void MvtxRawEvtHeader::AddFeeId( const std::set<int>& mvtxFeeIds)
{
  m_MvtxFeeIdSet.insert(mvtxFeeIds.cbegin(), mvtxFeeIds.cend());
}

void MvtxRawEvtHeader::AddL1Trg( const std::set<uint64_t>& mvtxL1TrgSet)
{
  m_MvtxL1TrgSet.insert(mvtxL1TrgSet.cbegin(), mvtxL1TrgSet.cend());
}
