#include "MbdPacketv1.h"

void MbdPacketv1::Reset()
{
  OfflinePacketv1::Reset();
  bunchnumber = std::numeric_limits<char>::min();
  return;
}
