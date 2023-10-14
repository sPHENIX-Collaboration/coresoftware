#include "MvtxRawEvtHeader.h"

void MvtxRawEvtHeader::Reset()
{
  m_MvtxFeeIdSet.clear();
}

void MvtxRawEvtHeader::identify(std::ostream &os) const
{
  os << "MvtxRawEvtHeader" << std::endl;
  os << "Feeids in the events:  " << m_MvtxFeeIdSet.size()<< std::endl;
}

int MvtxRawEvtHeader::isValid() const
{
  return  m_MvtxFeeIdSet.size();
}

void MvtxRawEvtHeader::AddFeeId( const std::set<int>& mvtxFeeIds)
{
  m_MvtxFeeIdSet.insert(mvtxFeeIds.cbegin(), mvtxFeeIds.cend());
}
