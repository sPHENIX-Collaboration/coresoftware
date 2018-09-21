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
#include <string>
#include <vector>
class PHCompositeNode;
class Fun4AllHistoManager;
class TTree;

class TPCFEETestRecov1 : public SubsysReco
{
 public:
  TPCFEETestRecov1(const std::string &outputfilename =
                       "TPCFEETestRecov1.root");
  virtual ~TPCFEETestRecov1();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

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

    uint32_t pedestal;
    uint32_t max;

    ChannelHeader()
      : size(0)
      , packet_type(0)
      , bx_counter(0)
      , sampa_address(0)
      , sampa_channel(0)
      , fee_channel(0)
      , pad_x(-1)
      , pad_y(-1)
      , pedestal(0)
      , max(0)
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

  TTree *m_chanT;

  ChannelHeader m_chanHeader;
  ChannelHeader *m_pchanHeader;  //! ->m_chanHeader,  for filling TTree
  std::vector<uint32_t> m_chanData;

  // clustering stuff

  //! full event data
  std::vector<std::vector<std::vector<int>>> m_data;

  void RoughZeroSuppression(std::vector<int> &data);
  void Clustering(void);

#endif  // #ifndef __CINT__
};

#endif /* CORESOFTWARE_OFFLINE_PACKAGES_TPCDAQ_TPCFEETESTRECOV1_H_ */
