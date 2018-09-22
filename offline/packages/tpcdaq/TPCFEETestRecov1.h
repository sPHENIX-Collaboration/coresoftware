/*
 * TPCFEETestRecov1.h
 *
 *  Created on: Sep 19, 2018
 *      Author: jinhuang
 */

#ifndef CORESOFTWARE_OFFLINE_PACKAGES_TPCDAQ_TPCFEETESTRECOV1_H_
#define CORESOFTWARE_OFFLINE_PACKAGES_TPCDAQ_TPCFEETESTRECOV1_H_

#include <fun4all/SubsysReco.h>

#include <TObject.h>

#include <stdint.h>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

class PHCompositeNode;
class Fun4AllHistoManager;
class TTree;
class TClonesArray;

namespace TPCDaqDefs
{
namespace FEEv1
{
class SampleFit_PowerLawDoubleExp_PDFMaker;
}
}  // namespace TPCDaqDefs

class TPCFEETestRecov1 : public SubsysReco
{
 public:
  TPCFEETestRecov1(const std::string &outputfilename =
                       "TPCFEETestRecov1.root");
  virtual ~TPCFEETestRecov1();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int ResetEvent(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void setClusteringZeroSuppression(int threshold)
  {
    m_clusteringZeroSuppression = threshold;
  }

  void setNPostSample(int nPostSample)
  {
    m_nPostSample = nPostSample;
  }

  void setNPreSample(int nPreSample)
  {
    m_nPreSample = nPreSample;
  }

  //! simple event header class for ROOT file IO
  class EventHeader : public TObject
  {
   public:
    int run;
    int event;

    uint32_t bx_counter;

    EventHeader()
      : run(-1)
      , event(-1)
      , bx_counter(0)
    {
    }

    ClassDef(TPCFEETestRecov1::EventHeader, 1)
  };

  //! buffer for full event data
  class PadPlaneData
  {
   public:
    PadPlaneData();
    void Reset();

    struct SampleID
    {
      int pady;
      int padx;
      int sample;

      void adjust(const SampleID &adjustment)
      {
        pady += adjustment.pady;
        padx += adjustment.padx;
        sample += adjustment.sample;
      }
    };

    static bool IsValidPad(const int pad_x, const int pad_y);
    std::vector<int> &getPad(const int pad_x, const int pad_y);
    int getSample(const SampleID &id);

    //! 3-D Graph clustering based on PHMakeGroups()
    void Clustering(int zero_suppression, bool verbosity = false);

#ifndef __CINT__

    const std::vector<std::vector<std::vector<int>>> &getData() const
    {
      return m_data;
    }

    const std::multimap<int, SampleID> &getGroups() const
    {
      return m_groups;
    }

   private:
    //! full event data in index order of m_data[pady][padx][sample]
    std::vector<std::vector<std::vector<int>>> m_data;

    std::multimap<int, SampleID> m_groups;

#endif  // #ifndef __CINT__
  };

  //! buffer for a cluster's data
  class ClusterData : public TObject
  {
   public:
    ClusterData()
      : min_sample(-1)
      , max_sample(-1)
      , peak(NAN)
      , peak_sample(NAN)
      , pedstal(NAN)
      , avg_padx(NAN)
      , avg_pady(NAN)
    {
    }

    std::set<int> padxs;
    std::set<int> padys;
    std::set<int> samples;

#ifndef __CINT__
    std::map<int, std::vector<double>> padx_samples;
    std::map<int, std::vector<double>> pady_samples;
    std::vector<double> sum_samples;
#endif  // #ifndef __CINT__

    int min_sample;
    int max_sample;

    double peak;
    double peak_sample;
    double pedstal;

    std::map<int, double> padx_peaks;
    std::map<int, double> pady_peaks;

    double avg_padx;
    double avg_pady;

    ClassDef(ClusterData, 1);
  };

  //! simple channel header class for ROOT file IO
  class ChannelHeader : public TObject
  {
   public:
    int size;
    //! = p->iValue(channel * kPACKET_LENGTH + 2) & 0xffff;  // that's the Elink packet type
    uint8_t packet_type;
    //! = ((p->iValue(channel * kPACKET_LENGTH + 4) & 0xffff) << 4) | (p->iValue(channel * kPACKET_LENGTH + 5) & 0xffff);
    uint32_t bx_counter;
    //! = (p->iValue(channel * kPACKET_LENGTH + 3) >> 5) & 0xf;
    uint8_t sampa_address;
    //! = p->iValue(channel * kPACKET_LENGTH + 3) & 0x1f;
    uint16_t sampa_channel;
    //! = (sampa_address << 5) | sampa_channel;
    uint16_t fee_channel;

    //! pad coordinate
    int pad_x;
    int pad_y;

    int pedestal;
    int max;

    ChannelHeader()
      : size(0)
      , packet_type(0)
      , bx_counter(0)
      , sampa_address(0)
      , sampa_channel(0)
      , fee_channel(0)
      , pad_x(-1)
      , pad_y(-1)
      , pedestal(-1)
      , max(-1)
    {
    }

    ClassDef(TPCFEETestRecov1::ChannelHeader, 1)
  };

 private:
#ifndef __CINT__

  // IO stuff

  Fun4AllHistoManager *getHistoManager();

  std::string m_outputFileName;

  TTree *m_eventT;

  EventHeader m_eventHeader;
  EventHeader *m_peventHeader;  //! ->m_eventHeader,  for filling TTree

  int m_nClusters;
  TClonesArray *m_IOClusters;

  TTree *m_chanT;

  ChannelHeader m_chanHeader;
  ChannelHeader *m_pchanHeader;  //! ->m_chanHeader,  for filling TTree
  std::vector<uint32_t> m_chanData;

  // clustering stuff
  PadPlaneData m_padPlaneData;
  std::map<int, ClusterData> m_clusters;

  //! rough zero suppression by subtracting sample medium value
  //! \return pair of pedestal and max-pedestal
  static std::pair<int, int> roughZeroSuppression(std::vector<int> &data);

  //! Clustering then prepare IOs
  void Clustering(void);

#endif  // #ifndef __CINT__

  int m_clusteringZeroSuppression;
  int m_nPreSample;
  int m_nPostSample;
  TPCDaqDefs::FEEv1::SampleFit_PowerLawDoubleExp_PDFMaker *m_pdfMaker;
};

bool operator<(const TPCFEETestRecov1::PadPlaneData::SampleID &s1, const TPCFEETestRecov1::PadPlaneData::SampleID &s2);

#endif /* CORESOFTWARE_OFFLINE_PACKAGES_TPCDAQ_TPCFEETESTRECOV1_H_ */
