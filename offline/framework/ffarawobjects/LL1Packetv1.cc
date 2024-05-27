#include "LL1Packetv1.h"

#include <Event/packetConstants.h>

#include <TSystem.h>

#include <iomanip>

LL1Packetv1::LL1Packetv1()
{
  for (auto &row : samples)
  {
    row.fill(0);
  }
}

void LL1Packetv1::Reset()
{
  OfflinePacketv1::Reset();
  PacketEvtSequence = 0;
  NrChannels = 0;
  NrSamples = 0;
  TriggerWords = 0;
  SlotNr = 0;
  CardNr = 0;
  for (auto &row : samples)
  {
    row.fill(0);
  }
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
  if (what == "Channels")
  {
    return getNrChannels();
  }
  if (what == "TRIGGERWORDS")
  {
    return getTriggerWords();
  }
  if (what == "SLOTNR")
  {
    return getSlotNr();
  }
  if (what == "CARDNR")
  {
    return getCardNr();
  }
  if (what == "MONITOR")
  {
    return getMonitor();
  }

  std::cout << "invalid selection " << what << std::endl;
  return std::numeric_limits<int>::min();
}

int LL1Packetv1::iValue(const int channel, const int sample) const
{
  return samples.at(channel).at(sample);
}

void LL1Packetv1::identify(std::ostream &os) const
{
  os << "LL1Packetv1: " << std::endl;
  OfflinePacketv1::identify(os);
  os << "Pkt Event no: " << getPacketEvtSequence() << std::endl;
  os << "Hitformat: " << getHitFormat() << std::endl;
  os << "Number of Channels: " << getNrChannels() << std::endl;
  os << "Number of Samples: " << getNrSamples() << std::endl;
}

void LL1Packetv1::dump(std::ostream &os) const
{
  os << std::dec << std::setprecision(2) << "Trigger Module = " << iValue(0, "SLOTNR") * 2 + iValue(0, "CARDNR") << std::endl;
  os << std::dec << std::setprecision(4) << "Evt Nr = " << iValue(0, "EVTNR") << std::endl;
  os << std::dec << std::setprecision(4) << "Clock = " << iValue(0, "CLOCK") << std::endl;
  os << std::dec << std::setprecision(4) << "Monitor = " << iValue(0, "MONITOR") << std::endl;

  switch (getHitFormat())
  {
  case IDLL1_MBD:
    dump_idll1_mbd(os);
    break;
  case IDLL1_EMCAL_MON3:
    dump_idll_emcal_mon3(os);
    break;
  case IDLL1_JET_EMCAL_MON1:
    dump_idll_jet_emcal_mon1(os);
    break;
  default:
    std::cout << "unknown hit format: "
              << getHitFormat() << std::endl;
    //    gSystem->Exit(1);
  }
  return;
}

void LL1Packetv1::dump_idll1_mbd(std::ostream &os) const
{
  for (int ifem = 0; ifem < 4; ifem++)
  {
    std::cout << std::dec << "FEM : " << ifem << std::endl;
    for (int iq = 0; iq < 8; iq++)
    {
      std::cout << std::dec << "Q" << iq << "\t||  \t";
      for (int is = 0; is < iValue(0, "SAMPLES"); is++)
      {
        std::cout << std::hex << iValue(is, ifem * 13 + iq) << "\t";
      }
      std::cout << " |" << std::endl;
    }
    std::cout << std::dec << "NH \t||  \t";
    for (int is = 0; is < iValue(0, "SAMPLES"); is++)
    {
      std::cout << std::hex << iValue(is, ifem * 13 + 8) << "\t";
    }
    std::cout << " |" << std::endl;

    for (int iq = 0; iq < 4; iq++)
    {
      std::cout << std::dec << "T" << iq << "\t||  \t";
      for (int is = 0; is < iValue(0, "SAMPLES"); is++)
      {
        std::cout << std::hex << iValue(is, ifem * 13 + 9 + iq) << "\t";
      }
      std::cout << " |" << std::endl;
    }
    std::cout << " " << std::endl;
  }

  for (int iw = 0; iw < iValue(0, "TRIGGERWORDS"); iw++)
  {
    std::cout << std::dec << "W " << iw << "\t||  \t";
    for (int is = 0; is < iValue(0, "SAMPLES"); is++)
    {
      std::cout << std::hex << iValue(is, 52 + iw) << "\t";
    }

    std::cout << " |" << std::endl;
  }
}
//  Packet_iddigitizerv3
void LL1Packetv1::dump_idll_emcal_mon3(std::ostream &os) const
{
  std::cout << "-------------------------------------------------------------- " << std::endl;
  for (int ch = 0; ch < 24; ch++)
  {
    std::cout << std::dec << "Fiber: " << ch << std::endl;
    for (int ic = 0; ic < 16; ic++)
    {
      std::cout << std::dec << " Sum " << ic << " |";
      for (int is = 0; is < iValue(0, "SAMPLES"); is++)
      {
        std::cout << std::hex << " " << iValue(is, ch * 16 + ic);
      }

      std::cout << " |" << std::endl;
      std::cout << "-------------------------------------------------------------- " << std::endl;
    }
  }
  for (int is = 0; is < iValue(0, "SAMPLES"); is++)
  {
    std::cout << std::dec << "Sample: " << is << std::endl;
    for (int ic = 0; ic < 2; ic++)
    {
      for (int ie = 0; ie < 12; ie++)
      {
        std::cout << std::hex << " " << iValue(is, 24 * 16 + ic * 12 + ie);
      }

      std::cout << " |" << std::endl;
      std::cout << "-------------------------------------------------------------- " << std::endl;
    }
  }
}
void LL1Packetv1::dump_idll_jet_emcal_mon1(std::ostream &os) const
{
  std::cout << " -------------- " << (iValue(0, "MONITOR") ? "HCAL Data Map" : "EMCAL Data Map") << " -------------- " << std::endl;
  for (int sample = 0; sample < iValue(0, "SAMPLES"); sample++)
  {
    std::cout << std::dec << "BC : " << sample << std::endl;
    std::cout << std::dec << "\t phibin --> " << std::endl;
    std::cout << std::dec << "etabin\t||  \t" << std::endl;
    for (int ic = 0; ic < 32; ic++) std::cout << std::dec << "\t" << ic;
    std::cout << " " << std::endl;
    for (int ic = 0; ic < 12; ic++)
    {
      std::cout << std::dec << ic << "\t||";
      for (int is = 0; is < 32; is++)
      {
        std::cout << std::hex << "\t" << iValue(sample, ic * 32 + is);
      }
      std::cout << " |" << std::endl;
    }
    std::cout << " " << std::endl;
  }

  for (int is = 0; is < iValue(0, "SAMPLES"); is++)
  {
    std::cout << std::dec << "Sample: " << is << std::endl;
    for (int ic = 0; ic < 9; ic++)
    {
      for (int ie = 0; ie < 32; ie++)
      {
        std::cout << std::hex << " " << iValue(is, 12 * 32 + ic * 32 + ie);
      }

      std::cout << " |" << std::endl;
      std::cout << "-------------------------------------------------------------- " << std::endl;
    }
  }
}
