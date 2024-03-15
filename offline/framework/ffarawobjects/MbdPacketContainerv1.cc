#include "MbdPacketContainerv1.h"
#include "MbdPacketv1.h"

#include <TClonesArray.h>

static const int NMBDPACKETS = 2;

MbdPacketContainerv1::MbdPacketContainerv1()
{
  MbdPacketsTCArray = new TClonesArray("MbdPacketv1", NMBDPACKETS);
}

MbdPacketContainerv1::~MbdPacketContainerv1()
{
  delete MbdPacketsTCArray;
}

void MbdPacketContainerv1::Reset()
{
  MbdPacketsTCArray->Clear();
  MbdPacketsTCArray->Expand(NMBDPACKETS);
}

void MbdPacketContainerv1::identify(std::ostream &os) const
{
  os << "MbdPacketContainerv1" << std::endl;
  os << "containing " << MbdPacketsTCArray->GetEntriesFast() << " Mbd hits" << std::endl;
  MbdPacket *mbdhit = static_cast<MbdPacket *>(MbdPacketsTCArray->At(0));
  if (mbdhit)
  {
    os << "for beam clock: " << std::hex << mbdhit->getBCO() << std::dec << std::endl;
  }
}

int MbdPacketContainerv1::isValid() const
{
  return MbdPacketsTCArray->GetSize();
}

unsigned int MbdPacketContainerv1::get_npackets()
{
  return MbdPacketsTCArray->GetEntriesFast();
}

MbdPacket *MbdPacketContainerv1::AddPacket()
{
  MbdPacket *newhit = new ((*MbdPacketsTCArray)[MbdPacketsTCArray->GetLast() + 1]) MbdPacketv1();
  return newhit;
}

MbdPacket *MbdPacketContainerv1::AddPacket(MbdPacket *mbdhit)
{
  // need a dynamic cast here to use the default copy ctor for MbdPacketv1
  // which copies the std::arrays
  MbdPacket *newhit = new ((*MbdPacketsTCArray)[MbdPacketsTCArray->GetLast() + 1]) MbdPacketv1(*(dynamic_cast<MbdPacketv1 *>(mbdhit)));
  return newhit;
}

MbdPacket *MbdPacketContainerv1::getPacket(unsigned int index)
{
  return (MbdPacket *) MbdPacketsTCArray->At(index);
}
