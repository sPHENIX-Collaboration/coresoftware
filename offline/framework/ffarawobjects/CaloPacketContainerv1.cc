#include "CaloPacketContainerv1.h"
#include "CaloPacketv1.h"

#include <TClonesArray.h>

static const int NCALOPACKETS = 2;

CaloPacketContainerv1::CaloPacketContainerv1()
{
  CaloPacketsTCArray = new TClonesArray("CaloPacketv1",NCALOPACKETS);
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
  os << "containing " << CaloPacketsTCArray->GetEntriesFast() << " Calo hits" << std::endl;
  CaloPacket *calohit = static_cast< CaloPacket *> (CaloPacketsTCArray->At(0));
  if (calohit)
   {
     os << "for beam clock: " << std::hex << calohit->getBCO() << std::dec << std::endl;
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
   CaloPacket *newhit = new ((*CaloPacketsTCArray)[CaloPacketsTCArray->GetLast()+1]) CaloPacketv1();
   return newhit;
 }

CaloPacket *CaloPacketContainerv1::AddPacket(CaloPacket *calohit)
 {
// need a dynamic cast here to use the default copy ctor for CaloPacketv1
// which copies the std::arrays
   CaloPacket *newhit = new ((*CaloPacketsTCArray)[CaloPacketsTCArray->GetLast()+1]) CaloPacketv1(*(dynamic_cast<CaloPacketv1 *> (calohit)));
   return newhit;
 }

CaloPacket *CaloPacketContainerv1::getPacket(unsigned int index)
{
   return (CaloPacket *) CaloPacketsTCArray->At(index);
}
