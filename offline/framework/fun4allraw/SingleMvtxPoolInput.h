#ifndef FUN4ALLRAW_SINGLEMVTXPOOLINPUT_H
#define FUN4ALLRAW_SINGLEMVTXPOOLINPUT_H

#include "SingleStreamingInput.h"

#include <algorithm>
#include <map>
#include <vector>

class MvtxRawHit;
class Packet;
class mvtx_pool;

class SingleMvtxPoolInput : public SingleStreamingInput
{
 public:
  explicit SingleMvtxPoolInput(const std::string &name);
  ~SingleMvtxPoolInput() override;
  void FillPool(const uint64_t minBCO) override;
  void CleanupUsedPackets(const uint64_t bclk) override;
  bool CheckPoolDepth(const uint64_t bclk) override;
  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;

  void SetBcoRange(const unsigned int i) { m_BcoRange = i; }
  void ConfigureStreamingInputManager() override;
  void SetNegativeBco(const unsigned int value) { m_NegativeBco = value; }
  void setRawEventHeaderName(const std::string &name) { m_rawEventHeaderName = name; }
  const std::string &getRawEventHeaderName() const { return m_rawEventHeaderName; }

  void  SetReadStrWidthFromDB(const bool val){ m_readStrWidthFromDB = val; }
  bool  GetReadStrWidthFromDB(){ return m_readStrWidthFromDB; }
  void  SetStrobeWidth(const float val) { m_strobeWidth = val; }
  float GetStrobeWidth() { return m_strobeWidth; }

 protected:
 private:
  Packet **plist{nullptr};
  unsigned int m_NumSpecialEvents{0};
  unsigned int m_BcoRange{0};
  unsigned int m_NegativeBco{0};
  std::string m_rawEventHeaderName = "MVTXRAWEVTHEADER";

  std::map<uint64_t, std::vector<MvtxRawHit *>> m_MvtxRawHitMap;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::map<int, uint64_t> m_FeeStrobeMap;
  std::set<uint64_t> m_BclkStack;
  std::set<uint64_t> gtmL1BcoSet;  // GTM L1 BCO
  std::map<int, mvtx_pool *> poolmap;

  bool m_readStrWidthFromDB = true;
  float m_strobeWidth = 0;
};

#endif
