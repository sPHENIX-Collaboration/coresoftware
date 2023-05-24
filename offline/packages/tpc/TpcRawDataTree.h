// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TpcRawDataTree_H
#define TpcRawDataTree_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;

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
  void RemoveAllPackets()
  {
    m_packets.clear();
  }

 protected:
  //! which packet to decode
  std::vector<int> m_packets;

 private:

  std::string m_fname;
  TFile * m_file = nullptr;
  TTree * m_SampleTree = nullptr;
  TTree * m_PacketTree = nullptr;

  int m_packet = 0;
  int m_frame = 0;
  int m_nWaveormInFrame = 0;
  int m_maxFEECount = 0;
  int m_nSamples = 0;
  int m_fee = 0;
  int m_sampaAddress = 0;
  int m_sampaChannel = 0;
  int m_Channel = 0;
  int m_BCO = 0;
  int m_checksum = 0;
  int m_checksumError = 0;
  std::vector<unsigned short> m_adcSamples;
};

#endif  // TpcRawDataTree_H
