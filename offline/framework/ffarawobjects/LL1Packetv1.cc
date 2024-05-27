#include "LL1Packetv1.h"

#include <TSystem.h>

#include <iomanip>

LL1Packetv1::LL1Packetv1()
{
}

void LL1Packetv1::Reset()
{
  OfflinePacketv1::Reset();
  PacketEvtSequence = 0;
  NrChannels = 0;
  NrSamples = 0;
  return;
}

int LL1Packetv1::iValue(const int n, const std::string &what) const
{
  if (what == "CLOCK")
  {
    return getBCO();
  }

  if (what == "EVTNR")
  {
    return getPacketEvtSequence();
  }

  if (what == "SAMPLES")
  {
    return getNrSamples();
  }

  std::cout << "invalid selection " << what << std::endl;
  return std::numeric_limits<int>::min();
}

 int LL1Packetv1::iValue(const int channel, const int sample) const
 {
   if (channel > 10)
   {
     return sample;
   }
   return 0;
//   return samples.at(channel).at(sample);
 }

void LL1Packetv1::identify(std::ostream &os) const
{
  os << "LL1Packetv1: " << std::endl;
  OfflinePacketv1::identify(os);
  os << "Pkt Event no: " << getPacketEvtSequence() << std::endl;
  os << "FEM Event no: " << std::hex;
}

void LL1Packetv1::dump(std::ostream &os) const
{
  switch (getHitFormat())
  {
  case 93:
    dump93(os);
    break;
  case 172:
    dump172(os);
    break;
  default:
    std::cout << "unknown hit format: "
	      << getHitFormat() << std::endl;
    gSystem->Exit(1);
  }
  return;
}


void LL1Packetv1::dump93(std::ostream &os) const
{
  int _nchannels = iValue(0, "CHANNELS");
  int _nsamples = iValue(0, "SAMPLES");
  os << "Evt Nr:      " << iValue(0, "EVTNR") << std::endl;
  os << "Clock:       " << iValue(0, "CLOCK") << std::endl;
  os << "Channels:    " << iValue(0, "CHANNELS") << std::endl;
  os << "Samples:     " << iValue(0, "SAMPLES") << std::endl;
}
//  Packet_iddigitizerv3
void LL1Packetv1::dump172(std::ostream &os) const
{
  int _nchannels = iValue(0, "CHANNELS");
  int _nsamples = iValue(0, "SAMPLES");
  os << "Evt Nr:      " << iValue(0, "EVTNR") << std::endl;
  os << "Clock:       " << iValue(0, "CLOCK") << std::endl;
  os << "Channels:    " << iValue(0, "CHANNELS") << std::endl;
  os << "Samples:     " << iValue(0, "SAMPLES") << std::endl;
}
