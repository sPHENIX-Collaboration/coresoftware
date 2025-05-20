#include "LL1PacketContainerv1.h"
#include "LL1Packetv1.h"

#include <phool/phool.h>

#include <TClonesArray.h>

static const int NLL1PACKETS = 18;

LL1PacketContainerv1::LL1PacketContainerv1()
  : LL1PacketsTCArray(new TClonesArray("LL1Packetv1", NLL1PACKETS))
{
}

LL1PacketContainerv1::~LL1PacketContainerv1()
{
  delete LL1PacketsTCArray;
}

void LL1PacketContainerv1::Reset()
{
  LL1PacketsTCArray->Clear();
  LL1PacketsTCArray->Expand(NLL1PACKETS);
}

void LL1PacketContainerv1::identify(std::ostream &os) const
{
  os << "LL1PacketContainerv1" << std::endl;
  os << "containing " << LL1PacketsTCArray->GetEntriesFast() << " LL1 Packets" << std::endl;
  for (int i = 0; i <= LL1PacketsTCArray->GetLast(); i++)
  {
    LL1Packet *ll1pkt = dynamic_cast<LL1Packet *>(LL1PacketsTCArray->At(i));
    if (ll1pkt)
    {
      os << "id: " << ll1pkt->getIdentifier() << std::endl;
      os << "for beam clock: " << std::hex << ll1pkt->getBCO() << std::dec << std::endl;
    }
  }
}

int LL1PacketContainerv1::isValid() const
{
  return LL1PacketsTCArray->GetSize();
}

unsigned int LL1PacketContainerv1::get_npackets()
{
  return LL1PacketsTCArray->GetEntriesFast();
}

LL1Packet *LL1PacketContainerv1::AddPacket()
{
  LL1Packet *newhit = new ((*LL1PacketsTCArray)[LL1PacketsTCArray->GetLast() + 1]) LL1Packetv1();
  return newhit;
}

LL1Packet *LL1PacketContainerv1::AddPacket(LL1Packet *ll1packet)
{
  // need a dynamic cast here to use the default copy ctor for LL1Packetv1
  // which copies the std::arrays
  LL1Packet *newhit = new ((*LL1PacketsTCArray)[LL1PacketsTCArray->GetLast() + 1]) LL1Packetv1(*(dynamic_cast<LL1Packetv1 *>(ll1packet)));
  return newhit;
}

LL1Packet *LL1PacketContainerv1::getPacket(unsigned int index)
{
  return (LL1Packet *) LL1PacketsTCArray->At(index);
}

LL1Packet *LL1PacketContainerv1::getPacketbyId(int id)
{
  for (int i = 0; i <= LL1PacketsTCArray->GetLast(); i++)
  {
    LL1Packet *pkt = (LL1Packet *) LL1PacketsTCArray->At(i);
    if (pkt->getIdentifier() == id)
    {
      return pkt;
    }
  }
  return nullptr;
}
