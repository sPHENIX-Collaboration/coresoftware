// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPC_TPCRAWDATATREE_H
#define TPC_TPCRAWDATATREE_H

#include <tpc/TpcMap.h>

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;
class TH1;
class TH2;

//! Dump TPC raw data in PRDF format to a TTree for online debugging and seeding formal Fun4All reco/calib modules
class TpcRawDataTree : public SubsysReco
{
 public:
  explicit TpcRawDataTree(const std::string &fname = "TpcRawDataTree.root");

  ~TpcRawDataTree() override {}

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

  void AddPacket(int packet)
  {
    m_packets.push_back(packet);
  }
  void includeXYPos(bool doInclude)
  {
    m_includeXYPos = doInclude;
  }

 private:
  TH1 *R1_hist{nullptr};
  TH1 *R2_hist{nullptr};
  TH1 *R3_hist{nullptr};
  TH1 *TotalFEE{nullptr};
  TH1 *TotalFEEsampa{nullptr};
  TH1 *TotalFRAME{nullptr};
  TH1 *checksumError_fee{nullptr};
  TH1 *checksumError_feesampa{nullptr};
  TH1 *checksumError_frame{nullptr};

  TH2 *R1_time{nullptr};
  TH2 *R2_time{nullptr};
  TH2 *R3_time{nullptr};

  TFile *m_file{nullptr};

  TTree *m_SampleTree{nullptr};
  TTree *m_PacketTree{nullptr};
  TTree *m_TaggerTree{nullptr};

  uint64_t m_bco{0};
  uint64_t m_last_bco{0};

  double m_xPos{0.};
  double m_yPos{0.};

  uint32_t m_lvl1_count{0};
  uint32_t m_endat_count{0};

  int m_packet{0};
  int m_frame{0};
  int m_nWaveormInFrame{0};
  int m_maxFEECount{0};
  int m_nSamples{0};
  int m_fee{0};
  int m_sampaAddress{0};
  int m_sampaChannel{0};
  int m_Channel{0};
  int m_BCO{0};
  int m_checksum{0};
  int m_checksumError{0};
  int side{0};
  int m_nTaggerInFrame{0};
  int FEE_R[26]{2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3, 2, 2, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3};
  int FEE_map[26]{3, 2, 5, 3, 4, 0, 2, 1, 3, 4, 5, 7, 6, 2, 0, 1, 0, 1, 4, 5, 11, 9, 10, 8, 6, 7};

  uint16_t m_tagger_type{0};

  uint8_t m_is_endat{0};
  uint8_t m_is_lvl1{0};
  uint8_t m_modebits{0};

  bool m_includeXYPos{true};

  std::string m_fname;
  std::string sectorNum;

  //! which packet to decode
  std::vector<int> m_packets;
  std::vector<unsigned short> m_adcSamples;

  TpcMap M;
};

#endif  // TpcRawDataTree_H
