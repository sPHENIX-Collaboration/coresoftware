#include "OfflinePacketv1.h"

void OfflinePacketv1::Reset()
{
  evtseq = std::numeric_limits<int>::min();
  packetid = std::numeric_limits<int>::min();
  bco = std::numeric_limits<uint64_t>::max();
}
