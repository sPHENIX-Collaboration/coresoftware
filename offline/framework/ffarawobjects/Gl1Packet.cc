#include "Gl1Packet.h"

#include <phool/phool.h>

#include <iostream>

int Gl1Packet::iValue(const int i) const
{
  if (i == 0)
  {
    return getPacketNumber();
  }
  std::cout << PHWHERE << " Bad argument for iValue: " << i << std::endl;
  return std::numeric_limits<int>::min();
}

long long Gl1Packet::lValue(const int i, const int j) const
{
  return getScaler(i, j);
}

long long Gl1Packet::lValue(const int i, const std::string &what) const
{
  if (what == "BCO")
  {
    return getBCO();
  }
  if (what == "TriggerInput")
  {
    return getTriggerInput();
  }
  if (what == "TriggerVector")
  {
    return getTriggerVector();
  }
  if (what == "LiveVector")
  {
    return getLiveVector();
  }
  if (what == "ScaledVector")
  {
    return getScaledVector();
  }
  if (what == "GTMBusyVector")
  {
    return getGTMBusyVector();
  }
  if (what == "GTMAllBusyVector")
  {
    return getGTMAllBusyVector();
  }
  if (what == "BunchNumber")
  {
    return getBunchNumber();
  }
  if (what == "GL1PRAW")
  {
    return getGl1pScaler(i, 0);
  }
  if (what == "GL1PLIVE")
  {
    return getGl1pScaler(i, 1);
  }
  if (what == "GL1PSCALED")
  {
    return getGl1pScaler(i, 2);
  }
  std::cout << "option " << what << " not implemented" << std::endl;
  return std::numeric_limits<uint64_t>::max();
}
