#ifndef BCOLUMICOUNT_STREAMINGLUMIRECO_H
#define BCOLUMICOUNT_STREAMINGLUMIRECO_H

#include "StreamingLumiInfo.h"

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllHistoManager.h>

#include <string>
#include <utility>
#include <array>

class TH1;

class StreamingLumiReco : public SubsysReco
{
 public:
  StreamingLumiReco(const std::string &name = "STREAMINGBCOLUMIRECO");
  ~StreamingLumiReco() override = default;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;

  virtual const std::array<double, 120> get_bunchnumber_lumi() const { return m_bunchnumber_lumi; }

  virtual double get_lumi() const { return m_lumi; }


 private:
  static int CreateNodeTree(PHCompositeNode *topNode);

  int m_bunches = 120;

  //double m_xsec_MBDNS = 24.07*1e9; //convert to pb from Vernier scan DOUBLE CHECK VALUE!

  uint64_t m_rawgl1scaler{0};

  std::array<long, 120> m_bunchnumber_MBDNS{0};
  //std::array<long, 120> m_bunchnumber_MBDNS_live{0};
  //std::array<long, 120> m_bunchnumber_MBDNS_scaled{0};

  std::array<long, 120> m_bunchnumber_crossings{0};

  std::array<double, 120> m_bunchnumber_lumi{0.};

  double m_lumi{0.};

  StreamingLumiInfo *m_streaming_lumi_info = nullptr;

};

#endif  // BCOLUMICOUNT_STREAMINGLUMIRECO_H
