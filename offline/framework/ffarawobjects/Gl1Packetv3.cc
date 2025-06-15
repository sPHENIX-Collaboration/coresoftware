#include "Gl1Packetv3.h"

#include <phool/phool.h>

#include <iomanip>

void Gl1Packetv3::Reset()
{
  Gl1Packetv2::Reset();
  GTMAllBusyVector = 0;
  return;
}

void Gl1Packetv3::identify(std::ostream &os) const
{
  os << "Gl1Packetv3: " << std::endl;
  os << "Id: " << getIdentifier() << std::endl;
  os << "EvtSeq: " << getEvtSequence() << std::endl;
  os << "BCO: 0x" << std::hex << getBCO() << std::dec << std::endl;
  os << "bunch number: " << getBunchNumber() << std::endl;
  return;
}

void Gl1Packetv3::FillFrom(const Gl1Packet *pkt)
{
  Gl1Packetv2::FillFrom(pkt);
  setGTMAllBusyVector(pkt->getGTMAllBusyVector());
}

void Gl1Packetv3::dump(std::ostream &os) const
{
  os << "packet nr:       " << iValue(0) << std::endl;
  os << "Beam Clock:      "
     << "0x" << std::hex << Gl1Packetv2::lValue(0, "BCO") << std::dec << "   " << lValue(0, "BCO") << std::endl;
  os << "Trigger Input:   "
     << "0x" << std::hex << lValue(0, "TriggerInput") << std::dec << "   " << lValue(0, "TriggerInput") << std::endl;
  os << "Live Vector:     "
     << "0x" << std::hex << lValue(0, "LiveVector") << std::dec << "   " << lValue(0, "LiveVector") << std::endl;
  os << "Scaled Vector:   "
     << "0x" << std::hex << lValue(0, "ScaledVector") << std::dec << "   " << lValue(0, "ScaledVector") << std::endl;
  os << "GTM Busy Vector: "
     << "0x" << std::hex << lValue(0, "GTMBusyVector") << std::dec << "   " << lValue(0, "GTMBusyVector") << std::endl;
  os << "GTM All Busy Vector: "
     << "0x" << std::hex << lValue(0, "GTMAllBusyVector") << std::dec << "   " << lValue(0, "GTMAllBusyVector") << std::endl;
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
