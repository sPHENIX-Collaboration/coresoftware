#ifndef FUN4ALLRAW_SINGLEMICROMEGASPOOLINPUT_H
#define FUN4ALLRAW_SINGLEMICROMEGASPOOLINPUT_H

#include "SingleStreamingInput.h"

#include <array>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

class Fun4AllEvtInputPoolManager;
class MicromegasRawHit;
class Packet;

class SingleMicromegasPoolInput : public SingleStreamingInput
{
 public:
  explicit SingleMicromegasPoolInput(const std::string &name = "SingleMicromegasPoolInput");
  ~SingleMicromegasPoolInput() override;
  void FillPool(const unsigned int nevents = 1) override;
  void CleanupUsedPackets(const uint64_t bclk) override;
  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;

  void SetBcoRange(const unsigned int value) { m_BcoRange = value; }
  void ConfigureStreamingInputManager() override;
  void SetNegativeBco(const unsigned int value) { m_NegativeBco = value; }

 private:
  Packet **plist{nullptr};
  unsigned int m_NumSpecialEvents{0};
  unsigned int m_BcoRange{0};
  unsigned int m_NegativeBco{0};

  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<uint64_t, std::vector<MicromegasRawHit *>> m_MicromegasRawHitMap;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::set<uint64_t> m_BclkStack;

  //! store relevant information for bco matching between lvl1 and fee.
  using m_bco_matching_pair_t = std::pair<unsigned int, uint64_t>;
  class bco_matching_information_t
  {
    public:

    //! first gtm bco (40 bits)
    /** it is needed to be able to convert gtm bco in a predicted fee bco */
    bool m_has_gtm_bco_first = false;
    uint64_t m_gtm_bco_first = 0;

    //! first fee bco (20 bits)
    /** it is needed to be able to convert gtm bco in a predicted fee bco */
    bool m_has_fee_bco_first = false;
    unsigned int m_fee_bco_first = 0;

    //! list of available bco
    std::list<uint64_t> m_gtm_bco_list;

    //! matching between fee bco and lvl1 bco
    std::list<m_bco_matching_pair_t> m_bco_matching_list;

    //! need to truncate bco matching list to some decent value
    void truncate( unsigned int /* maxsize */ );

    //! get predicted fee_bco from gtm_bco
    unsigned int get_predicted_fee_bco( uint64_t ) const;

  };

  //! map bco_information_t to packet id
  using bco_matching_information_map_t = std::map<unsigned int, bco_matching_information_t>;
  bco_matching_information_map_t m_bco_matching_information_map;

};

#endif
