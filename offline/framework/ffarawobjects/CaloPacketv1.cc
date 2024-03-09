#include "CaloPacketv1.h"

void CaloPacketv1::Reset()
{
  OfflinePacketv1::Reset();
  PacketEvtSequence = 0;
  femevt.fill(std::numeric_limits<uint32_t>::max());
  femclock.fill(std::numeric_limits<uint32_t>::max());
  for (auto &row :  samples)
  {
    row.fill(std::numeric_limits<uint32_t>::max());
  }
  return;
}

void CaloPacketv1::identify(std::ostream &os) const
{
  os << "CaloPacketv1: " << std::endl;
    OfflinePacketv1::identify(os);
  os << "Pkt Event no: " <<  getPacketEvtSequence() << std::endl;
  os << "FEM Event no: " << std::hex;
  for (const auto clk : femevt)
  { 
    std::cout << clk << " ";
  }
  std::cout << std::dec << std::endl;
  os << "FEM clk: " << std::hex;
  for (const auto clk : femclock)
  { 
    std::cout << clk << " ";
  }
  std::cout << std::dec << std::endl;
/*
  for (auto &iter :  samples)
  {
    for (auto &iter2 : iter)
    {
    std::cout << "sample: " << iter2 << std::endl;
    }
  }
*/
}
