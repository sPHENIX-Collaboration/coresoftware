#ifndef FUN4ALLRAW_TPCTIMEFRAMEBUILDERBASE_H
#define FUN4ALLRAW_TPCTIMEFRAMEBUILDERBASE_H

#include <cstdint>
#include <string>
#include <vector>

class Packet;
class TpcRawHit;

class TpcTimeFrameBuilderBase
{
 public:
  virtual ~TpcTimeFrameBuilderBase() = default;

  virtual int ProcessPacket(Packet *) = 0;
  virtual bool isMoreDataRequired(const uint64_t &gtm_bco) const = 0;
  virtual void CleanupUsedPackets(const uint64_t &bclk) = 0;
  virtual std::vector<TpcRawHit *> &getTimeFrame(const uint64_t &gtm_bco) = 0;

  virtual void setVerbosity(int i) = 0;
  virtual void fillBadFeeMap() = 0;
  virtual void SaveDigitalCurrentDebugTTree(const std::string &name) = 0;
  virtual void SaveBXCounterSyncCDBTTree(const std::string &name) = 0;
};

#endif
