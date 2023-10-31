#ifndef FUN4ALLRAW_SINGLEMICROMEGASPOOLINPUT_H
#define FUN4ALLRAW_SINGLEMICROMEGASPOOLINPUT_H

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

class SingleMicromegasPoolInput : public SingleStreamingInput
{
 public:
  explicit SingleMicromegasPoolInput(const std::string &name = "SingleMicromegasPoolInput" );
  ~SingleMicromegasPoolInput() override;
  void FillPool(const unsigned int nevents = 1) override;
  void CleanupUsedPackets(const uint64_t bclk) override;
  // bool CheckPoolDepth(const uint64_t bclk) override;
  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;

  void SetBcoRange(const unsigned int value ) {m_BcoRange = value;}
  void ConfigureStreamingInputManager() override;

 private:
  Packet **plist = nullptr;
  unsigned int m_NumSpecialEvents = 0;
  unsigned int m_BcoRange = 0;
  
  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<uint64_t, std::vector<MicromegasRawHit *>> m_MicromegasRawHitMap;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::set<uint64_t> m_BclkStack;
  
  //! keep track of matching between fee and gtm_bco 
  class bco_alignment_t
  {
    public:
    
    //! available gtm bcos
    std::list<uint64_t> gtm_bco_list;

    //! current fee bco
    unsigned int fee_bco = 0;
    
    //! current gtm bco
    uint64_t gtm_bco = 0;      
  };

  //! max number of fees per single micromegas input
  static constexpr unsigned short m_max_fee = 26;

  //! keep one bco alignment object per fee
  std::array<bco_alignment_t, m_max_fee> m_bco_alignment_list;
  
};

#endif
