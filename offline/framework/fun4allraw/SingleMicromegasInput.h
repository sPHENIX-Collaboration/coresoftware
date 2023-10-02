#ifndef FUN4ALLRAW_SINGLEMICROMEGASINPUT_H
#define FUN4ALLRAW_SINGLEMICROMEGASINPUT_H

#include "SingleStreamingInput.h"

#include <array>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

class Eventiterator;
class Fun4AllEvtInputPoolManager;
class MicromegasRawHit;
class Packet;

class SingleMicromegasInput : public SingleStreamingInput
{
 public:
  explicit SingleMicromegasInput(const std::string &name = "SingleMicromegasInput" );
  ~SingleMicromegasInput() override;
  void FillPool(const unsigned int nevents = 1) override;
  void CleanupUsedPackets(const uint64_t bclk) override;
  bool CheckPoolDepth(const uint64_t bclk) override;
  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;

 private:
  Packet **plist = nullptr;
  unsigned int m_NumSpecialEvents = 0;

  //! map bco to packet
  std::map<unsigned int, uint64_t> m_packet_bco;
  
  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<uint64_t, std::vector<MicromegasRawHit *>> m_MicromegasRawHitMap;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::set<uint64_t> m_BclkStack;
};

#endif
