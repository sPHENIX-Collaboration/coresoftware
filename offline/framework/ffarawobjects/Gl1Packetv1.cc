#include "Gl1Packetv1.h"

void Gl1Packetv1::Reset()
{
  OfflinePacketv1::Reset();
  bunchnumber = std::numeric_limits<char>::min();
  return;
}
