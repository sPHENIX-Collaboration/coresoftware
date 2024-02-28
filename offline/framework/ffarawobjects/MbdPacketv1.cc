#include "MbdPacketv1.h"

void MbdPacketv1::Reset()
{
  OfflinePacketv1::Reset();
  femclock.fill(std::numeric_limits<uint32_t>::max());
  return;
}
