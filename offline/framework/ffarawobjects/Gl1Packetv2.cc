#include "Gl1Packetv2.h"

#include <phool/phool.h>

#include <iomanip>

void Gl1Packetv2::Reset()
{
  OfflinePacketv1::Reset();
  packet_nr = 0;
  BunchNumber = std::numeric_limits<uint64_t>::max();
  TriggerInput = 0;
  LiveVector = 0;
  ScaledVector = 0;
  GTMBusyVector = 0;
  for (auto &row : scaler)
  {
    row.fill(0);
  }
  return;
}

void Gl1Packetv2::identify(std::ostream &os) const
{
  os << "Gl1Packetv2: " << std::endl;
  OfflinePacketv1::identify(os);
  os << "bunch number: " << BunchNumber << std::endl;
  return;
}

void Gl1Packetv2::FillFrom(const Gl1Packet *pkt)
{
  setBunchNumber(pkt->getBunchNumber());
  setPacketNumber(pkt->getPacketNumber());
  setTriggerInput(pkt->getTriggerInput());
  setLiveVector(pkt->getLiveVector());
  setScaledVector(pkt->getScaledVector());
  setGTMBusyVector(pkt->getGTMBusyVector());
  for (int i = 0; i < 64; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      setScaler(i, j, pkt->lValue(i, j));
    }
  }
  std::string gl1p_names[3]{"GL1PRAW", "GL1PLIVE", "GL1PSCALED"};
  for (int i = 0; i < 16; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      setGl1pScaler(i, j, pkt->lValue(i, gl1p_names[j]));
    }
  }
  OfflinePacketv1::FillFrom(pkt);
}

void Gl1Packetv2::dump(std::ostream &os) const
{
  os << "packet nr:       " << iValue(0) << std::endl;
  os << "Beam Clock:      "
     << "0x" << std::hex << lValue(0, "BCO") << std::dec << "   " << lValue(0, "BCO") << std::endl;
  os << "Trigger Input:   "
     << "0x" << std::hex << lValue(0, "TriggerInput") << std::dec << "   " << lValue(0, "TriggerInput") << std::endl;
  os << "Live Vector:     "
     << "0x" << std::hex << lValue(0, "LiveVector") << std::dec << "   " << lValue(0, "LiveVector") << std::endl;
  os << "Scaled Vector:   "
     << "0x" << std::hex << lValue(0, "ScaledVector") << std::dec << "   " << lValue(0, "ScaledVector") << std::endl;
  os << "GTM Busy Vector: "
     << "0x" << std::hex << lValue(0, "GTMBusyVector") << std::dec << "   " << lValue(0, "GTMBusyVector") << std::endl;
  os << "Bunch Number:    " << lValue(0, "BunchNumber") << std::endl
     << std::endl;
  os << "Trg #                  raw              live              scaled" << std::endl;
  os << "----------------------------------------------------------------" << std::endl;

  int i;

  for (i = 0; i < 64; i++)
  {
    if (lValue(i, 0) || lValue(i, 1) || lValue(i, 2))
    {
      os << std::setw(3) << i << "    ";
      os << " " << std::setw(18) << lValue(i, 0)
         << " " << std::setw(18) << lValue(i, 1)
         << " " << std::setw(18) << lValue(i, 2)
         << std::endl;
    }
  }
  os << std::endl;
  os << "Gl1P #                raw              live              scaled" << std::endl;
  os << "----------------------------------------------------------------" << std::endl;

  for (i = 0; i < 16; i++)
  {
    if (lValue(i, "GL1PRAW") || lValue(i, "GL1PLIVE") || lValue(i, "GL1PSCALED"))
    {
      os << std::setw(3) << i << "    ";
      os << " " << std::setw(18) << lValue(i, "GL1PRAW")
         << " " << std::setw(18) << lValue(i, "GL1PLIVE")
         << " " << std::setw(18) << lValue(i, "GL1PSCALED")
         << std::endl;
    }
  }
  os << std::endl;
}
