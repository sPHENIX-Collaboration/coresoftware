#include "CaloPacketv1.h"

#include <iomanip>

void CaloPacketv1::Reset()
{
  OfflinePacketv1::Reset();
  PacketEvtSequence = 0;
  NrChannels = 0;
  NrSamples = 0;
  NrModules = 0;
  event_checksum = 0;
  odd_checksum = 0;
  calc_event_checksum = 0;
  calc_odd_checksum = 0;
  module_address = 0;
  detid = 0;

  femclock.fill(std::numeric_limits<uint32_t>::max());
  femevt.fill(std::numeric_limits<uint32_t>::max());
  femslot.fill(std::numeric_limits<uint32_t>::max());

  for (auto &row : samples)
  {
    row.fill(std::numeric_limits<uint32_t>::max());
  }
  return;
}

int CaloPacketv1::iValue(const int n, const std::string &what) const
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

  if (what == "NRMODULES")
  {
    return getNrModules();
  }

  if (what == "CHANNELS")
  {
    return getNrChannels();
  }

  if (what == "DETID")
  {
    return getDetId();
  }

  if (what == "MODULEADDRESS")
  {
    return getModuleAddress();
  }

  if (what == "FEMSLOT")
  {
    if (n < 0 || n >= getNrModules())
    {
      return 0;
    }
    return femslot.at(n);
  }

  if (what == "FEMEVTNR")
  {
    if (n < 0 || n >= getNrModules())
    {
      return 0;
    }
    return femevt.at(n);
  }

  if (what == "FEMCLOCK")
  {
    if (n < 0 || n >= getNrModules())
    {
      return 0;
    }
    return femclock.at(n);
  }

  if (what == "EVENCHECKSUM")
  {
    return getEvenChecksum();
  }

  if (what == "ODDCHECKSUM")
  {
    return getOddChecksum();
  }

  if (what == "CALCEVENCHECKSUM")
  {
    return getCalcEvenChecksum();
  }

  if (what == "CALCODDCHECKSUM")
  {
    return getCalcOddChecksum();
  }

  if (what == "EVENCHECKSUMOK")
  {
    if (getCalcEvenChecksum() < 0)
    {
      return -1;
    }
    if (getCalcEvenChecksum() == getEvenChecksum())
    {
      return 1;
    }
    return 0;
  }

  if (what == "ODDCHECKSUMOK")
  {
    if (getCalcOddChecksum() < 0)
    {
      return -1;
    }
    if (getCalcOddChecksum() == getOddChecksum())
    {
      return 1;
    }
    return 0;
  }

  if (what == "CHECKSUMOK")
  {
    if (getCalcOddChecksum() < 0 || getCalcEvenChecksum())
    {
      return -1;
    }
    if (getCalcEvenChecksum() == getEvenChecksum() &&
        getCalcOddChecksum() == getOddChecksum())
    {
      return 1;
    }
    return 0;
  }

  std::cout << "invalid selection " << what << std::endl;
  return std::numeric_limits<int>::min();
}

int CaloPacketv1::iValue(const int channel, const int sample) const
{
  return samples.at(channel).at(sample);
}

void CaloPacketv1::identify(std::ostream &os) const
{
  os << "CaloPacketv1: " << std::endl;
  OfflinePacketv1::identify(os);
  os << "Pkt Event no: " << getPacketEvtSequence() << std::endl;
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

void CaloPacketv1::dump(std::ostream &os) const
{
  int _nchannels = iValue(0, "CHANNELS");
  int _nsamples = iValue(0, "SAMPLES");
  os << "Evt Nr:      " << iValue(0, "EVTNR") << std::endl;
  os << "Clock:       " << iValue(0, "CLOCK") << std::endl;
  os << "Nr Modules:  " << iValue(0, "NRMODULES") << std::endl;
  os << "Channels:    " << iValue(0, "CHANNELS") << std::endl;
  os << "Samples:     " << iValue(0, "SAMPLES") << std::endl;
  os << "Det. ID:     " << std::hex << "0x" << iValue(0, "DETID") << std::dec << std::endl;
  os << "Mod. Addr:   " << std::hex << "0x" << iValue(0, "MODULEADDRESS") << std::dec << std::endl;

  os << "FEM Slot:    ";
  for (int i = 0; i < iValue(0, "NRMODULES"); i++)
  {
    os << std::setw(8) << iValue(i, "FEMSLOT");
  }
  os << std::endl;

  os << "FEM Evt nr:  ";
  for (int i = 0; i < iValue(0, "NRMODULES"); i++)
  {
    os << std::setw(8) << iValue(i, "FEMEVTNR");
  }
  os << std::endl;

  os << "FEM Clock:   ";
  for (int i = 0; i < iValue(0, "NRMODULES"); i++)
  {
    os << std::setw(8) << iValue(i, "FEMCLOCK");
  }
  os << std::endl;

  os << "Even chksum: 0x" << std::hex << iValue(0, "EVENCHECKSUM") << "   calculated:  0x" << iValue(0, "CALCEVENCHECKSUM");
  if (iValue(0, "EVENCHECKSUMOK") == 1)
  {
    os << " ok" << std::endl;
  }
  else if (iValue(0, "EVENCHECKSUMOK") == 0)
  {
    os << " **wrong" << std::endl;
  }

  os << "Odd chksum:  0x" << std::hex << iValue(0, "ODDCHECKSUM") << "   calculated:  0x" << iValue(0, "CALCODDCHECKSUM");
  if (iValue(0, "ODDCHECKSUMOK") == 1)
  {
    os << " ok" << std::endl;
  }
  else if (iValue(0, "ODDCHECKSUMOK") == 0)
  {
    os << " **wrong" << std::endl;
  }
  os << std::dec << std::endl;

  for (int c = 0; c < _nchannels; c++)
  {
    os << std::setw(4) << c << " | ";

    os << std::hex;
    for (int s = 0; s < _nsamples; s++)
    {
      os << std::setw(6) << iValue(s, c);
    }
    os << std::dec << std::endl;
  }
}
