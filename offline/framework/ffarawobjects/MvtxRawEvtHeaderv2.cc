#include "MvtxRawEvtHeaderv2.h"

#include "MvtxFeeIdInfov1.h"

#include <TClonesArray.h>

static const uint8_t NMVTXFEEID = 48 * 3;

MvtxRawEvtHeaderv2::MvtxRawEvtHeaderv2()
  : m_MvtxFeeIdInfoTCArray(new TClonesArray("MvtxFeeIdInfov1", NMVTXFEEID))
{
}

MvtxRawEvtHeaderv2::~MvtxRawEvtHeaderv2()
{
  delete m_MvtxFeeIdInfoTCArray;
}

void MvtxRawEvtHeaderv2::Reset()
{
  m_MvtxL1TrgSet.clear();
  m_MvtxFeeIdInfoTCArray->Clear();
  m_MvtxFeeIdInfoTCArray->Expand(NMVTXFEEID);
}

void MvtxRawEvtHeaderv2::identify(std::ostream &os) const
{
  os << "MvtxRawEvtHeaderv2" << std::endl;
  os << "Feeids in the event:  " << m_MvtxFeeIdInfoTCArray->GetEntriesFast() << std::endl;
  os << "Mvtx L1 triggers in the event:  " << m_MvtxL1TrgSet.size() << std::endl;
}

int MvtxRawEvtHeaderv2::isValid() const
{
  return (!m_MvtxL1TrgSet.empty() || m_MvtxFeeIdInfoTCArray->GetEntriesFast());
}

MvtxFeeIdInfo *MvtxRawEvtHeaderv2::AddFeeIdInfo()
{
  MvtxFeeIdInfo *newFeeIdInfo = new ((*m_MvtxFeeIdInfoTCArray)[m_MvtxFeeIdInfoTCArray->GetLast() + 1]) MvtxFeeIdInfov1();
  return newFeeIdInfo;
}

MvtxFeeIdInfo *MvtxRawEvtHeaderv2::AddFeeIdInfo(MvtxFeeIdInfo *feeIdInfo)
{
  MvtxFeeIdInfo *newFeeIdInfo = new ((*m_MvtxFeeIdInfoTCArray)[m_MvtxFeeIdInfoTCArray->GetLast() + 1]) MvtxFeeIdInfov1(feeIdInfo);
  return newFeeIdInfo;
}

uint64_t MvtxRawEvtHeaderv2::get_nFeeIdInfo()
{
  return m_MvtxFeeIdInfoTCArray->GetEntriesFast();
}

MvtxFeeIdInfo *MvtxRawEvtHeaderv2::get_feeIdInfo(unsigned int index)
{
  return (MvtxFeeIdInfo *) m_MvtxFeeIdInfoTCArray->At(index);
}

void MvtxRawEvtHeaderv2::AddL1Trg(const std::set<uint64_t> &mvtxL1TrgSet)
{
  m_MvtxL1TrgSet.insert(mvtxL1TrgSet.cbegin(), mvtxL1TrgSet.cend());
}
