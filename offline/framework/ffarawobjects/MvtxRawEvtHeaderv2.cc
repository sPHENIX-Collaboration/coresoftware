#include "MvtxRawEvtHeaderv2.h"

void MvtxRawEvtHeaderv2::Reset()
{
  m_MvtxDetField.clear();
  m_MvtxL1TrgSet.clear();
}

void MvtxRawEvtHeaderv2::identify(std::ostream& os) const
{
  os << "MvtxRawEvtHeaderv2" << std::endl;
  os << "Feeids in the event:  " << m_MvtxDetField.size() << std::endl;
  os << "Mvtx L1 triggers in the event:  " << m_MvtxL1TrgSet.size() << std::endl;
}

int MvtxRawEvtHeaderv2::isValid() const
{
  return m_MvtxDetField.size();
}

void MvtxRawEvtHeaderv2::AddDetField(const std::map<uint16_t, uint32_t>& mvtxDetField)
{
  m_MvtxDetField.insert(mvtxDetField.cbegin(), mvtxDetField.cend());
}

void MvtxRawEvtHeaderv2::AddL1Trg(const std::set<uint64_t>& mvtxL1TrgSet)
{
  m_MvtxL1TrgSet.insert(mvtxL1TrgSet.cbegin(), mvtxL1TrgSet.cend());
}
