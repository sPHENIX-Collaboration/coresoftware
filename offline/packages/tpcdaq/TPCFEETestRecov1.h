/*
 * TPCFEETestRecov1.h
 *
 *  Created on: Sep 19, 2018
 *      Author: jinhuang
 */

#ifndef CORESOFTWARE_OFFLINE_PACKAGES_TPCDAQ_TPCFEETESTRECOV1_H_
#define CORESOFTWARE_OFFLINE_PACKAGES_TPCDAQ_TPCFEETESTRECOV1_H_

#include <fun4all/SubsysReco.h>

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

 private:
#ifndef __CINT__

  Fun4AllHistoManager *getHistoManager();

  std::string m_outputFileName;

  TTree *m_eventT;

  TTree *m_chanT;

  int m_event;

  struct ChannelHeader
  {
    uint16_t m_size;           //! = p->iValue(channel * kPACKET_LENGTH + 1) & 0xffff;        // number of words until the next channel (header included). this is the real packet_length
    uint8_t m_packet_type;     //! = p->iValue(channel * kPACKET_LENGTH + 2) & 0xffff;  // that's the Elink packet type
    uint32_t m_bx_counter;     //! = ((p->iValue(channel * kPACKET_LENGTH + 4) & 0xffff) << 4) | (p->iValue(channel * kPACKET_LENGTH + 5) & 0xffff);
    uint8_t m_sampa_address;   //! = (p->iValue(channel * kPACKET_LENGTH + 3) >> 5) & 0xf;
    uint16_t m_sampa_channel;  //! = p->iValue(channel * kPACKET_LENGTH + 3) & 0x1f;
    uint16_t m_fee_channel;    //! = (sampa_address << 5) | sampa_channel;

    ChannelHeader()
      : m_size(0)
      , m_packet_type(0)
      , m_bx_counter(0)
      , m_sampa_address(0)
      , m_sampa_channel(0)
      , m_fee_channel(0)
    {
    }
  };
  ChannelHeader m_chanHeader;
  std::vector<uint32_t> m_chanData;

#endif  // #ifndef __CINT__
};

#endif /* CORESOFTWARE_OFFLINE_PACKAGES_TPCDAQ_TPCFEETESTRECOV1_H_ */
