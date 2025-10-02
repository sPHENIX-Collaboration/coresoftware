// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPC_TPCMASKEDCHANNELMAP_H
#define TPC_TPCMASKEDCHANNELMAP_H

#include <tpc/TpcMap.h>
#include <trackbase/TrkrDefs.h>

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>
#include <set>
#include <TH2.h>

class PHCompositeNode;
class TpcRawHitContainer;

//! Dump TPC raw data in PRDF format to a TTree for online debugging and seeding formal Fun4All reco/calib modules
class TpcMaskedChannelMap : public SubsysReco
{
 public:
  explicit TpcMaskedChannelMap(const std::string &fname = "");

  ~TpcMaskedChannelMap() override {}

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
  */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
  */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void setHitCuts(float low, float high)
  {
    m_deadChanHitCut = low;
    m_hotChanHitCut = high;
  }

 private:
  
  int FEE_R[26]{2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3, 2, 2, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3};
  int FEE_map[26]{3, 2, 5, 3, 4, 0, 2, 1, 3, 4, 5, 7, 6, 2, 0, 1, 0, 1, 4, 5, 11, 9, 10, 8, 6, 7};
  int mc_sectors[12]{5, 4, 3, 2, 1, 0, 11, 10, 9, 8, 7, 6};
  float nhit_sectors_fees_channels[24][26][256] = {{{0}}};

  std::string m_run{""};
  std::string m_hotFile{""};
  std::string m_deadFile{""};

  TpcMap M;

  std::vector<TpcRawHitContainer*> rawhitcont_vec{};

  float m_deadChanHitCut = 0.1;
  float m_hotChanHitCut = 10;   
 
  float n_Events = 0; 

  std::string m_histogramFile{""};

  TH2F *h_hits_side0{nullptr};
  TH2F *h_hits_side1{nullptr};
  TH2F *h_masked_side0{nullptr};
  TH2F *h_masked_side1{nullptr};

  class channel_id_t
  {
    public:

    /// constructor
    channel_id_t( int layer, int sector, int side, int pad, float x, float y ):
      m_layer(layer),
      m_sector(sector),
      m_side(side),
      m_pad(pad),
      m_x(x),
      m_y(y)
    {}

    int m_layer = 0;
    int m_sector = 0;
    int m_side = 0;
    int m_pad = 0;
    int m_x = 0;
    int m_y = 0;

    // less than operator
    bool operator < (const channel_id_t& other ) const
    {
      if( m_layer != other.m_layer ) return m_layer < other.m_layer;
      else if( m_sector != other.m_sector ) return m_sector < other.m_sector;
      else if( m_side != other.m_side ) return m_side < other.m_side;
      else if( m_pad != other.m_pad ) return m_pad < other.m_pad;
      else if( m_x != other.m_x ) return m_x < other.m_x;
      else return m_y < other.m_y;
    }

    /// streamer
    friend std::ostream& operator << ( std::ostream& out, const channel_id_t& channel_id )
    {
      out << "{ " << channel_id.m_layer << ", " << channel_id.m_sector << ", " << channel_id.m_side << ", " << channel_id.m_pad << ", " << channel_id.m_x << ", " << channel_id.m_y << " }";
      return out;
    }

  };

  using channel_id_set_t = std::set<channel_id_t>;
  channel_id_set_t m_deadChannelCDB;
  channel_id_set_t m_hotChannelCDB;
};


#endif  // TpcMaskedChannelMap_H
