#include "CaloPacketContainerv1.h"
#include "CaloPacketv1.h"

#include <phool/phool.h>

#include <TClonesArray.h>

static const int NCALOPACKETS = 128;

CaloPacketContainerv1::CaloPacketContainerv1()
{
  CaloPacketsTCArray = new TClonesArray("CaloPacketv1", NCALOPACKETS);
}

CaloPacketContainerv1::~CaloPacketContainerv1()
{
  delete CaloPacketsTCArray;
}

void CaloPacketContainerv1::Reset()
{
  CaloPacketsTCArray->Clear();
  CaloPacketsTCArray->Expand(NCALOPACKETS);
}

void CaloPacketContainerv1::identify(std::ostream &os) const
{
  os << "CaloPacketContainerv1" << std::endl;
  os << "containing " << CaloPacketsTCArray->GetEntriesFast() << " Calo Packets" << std::endl;
  for (int i=0 ; i<=CaloPacketsTCArray->GetLast(); i++)
  {
    CaloPacket *calopkt = static_cast<CaloPacket *>(CaloPacketsTCArray->At(i));
    if (calopkt)
    {
      os << "id: " << calopkt->getIdentifier() << std::endl;
      os << "for beam clock: " << std::hex << calopkt->getBCO() << std::dec << std::endl;
    }
  }
}

int CaloPacketContainerv1::isValid() const
{
  return CaloPacketsTCArray->GetSize();
}

unsigned int CaloPacketContainerv1::get_npackets()
{
  return CaloPacketsTCArray->GetEntriesFast();
}

CaloPacket *CaloPacketContainerv1::AddPacket()
{
  CaloPacket *newhit = new ((*CaloPacketsTCArray)[CaloPacketsTCArray->GetLast() + 1]) CaloPacketv1();
  return newhit;
}

CaloPacket *CaloPacketContainerv1::AddPacket(CaloPacket *calohit)
{
  // need a dynamic cast here to use the default copy ctor for CaloPacketv1
  // which copies the std::arrays
  CaloPacket *newhit = new ((*CaloPacketsTCArray)[CaloPacketsTCArray->GetLast() + 1]) CaloPacketv1(*(dynamic_cast<CaloPacketv1 *>(calohit)));
  return newhit;
}

CaloPacket *CaloPacketContainerv1::getPacket(unsigned int index)
{
  return (CaloPacket *) CaloPacketsTCArray->At(index);
}

CaloPacket *CaloPacketContainerv1::getPacketbyId(int id)
{
  for (int i=0 ; i<=CaloPacketsTCArray->GetLast(); i++)
  {
    CaloPacket *pkt = (CaloPacket *) CaloPacketsTCArray->At(i);
    if (pkt->getIdentifier() == id)
    {
      return pkt;
    }
  }
  return nullptr;
}
