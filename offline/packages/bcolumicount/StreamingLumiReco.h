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

  virtual const std::array<double, 120> get_bunchnumber_lumi_raw() const { return m_bunchnumber_lumi_raw; }
  virtual const std::array<double, 120> get_bunchnumber_lumi_live() const { return m_bunchnumber_lumi_live; }
  virtual const std::array<double, 120> get_bunchnumber_lumi_scaled() const { return m_bunchnumber_lumi_scaled; }

  virtual double get_lumi_raw() const { return m_lumi_raw; }
  virtual double get_lumi_live() const { return m_lumi_live; }
  virtual double get_lumi_scaled() const { return m_lumi_scaled; }

  virtual void set_default_positive_window_length(int val) { m_default_positive_window_length = val; }
  virtual void set_default_negative_window_length(int val) { m_default_negative_window_length = val; }


 private:
  static int CreateNodeTree(PHCompositeNode *topNode);

  int m_bunches = 120;
  uint64_t m_bco{0};
  bool m_usable_bco_tag = false;
  std::pair<uint64_t, uint64_t> m_bco_streaming_window;
  unsigned int m_default_positive_window_length{340};
  unsigned int m_default_negative_window_length{20};

  //double m_xsec_MBDNS = 24.07*1e9; //convert to pb from Vernier scan DOUBLE CHECK VALUE!

  uint64_t m_rawgl1scaler{0};

  std::array<long, 120> m_bunchnumber_MBDNS_raw{0};
  std::array<long, 120> m_bunchnumber_MBDNS_live{0};
  std::array<long, 120> m_bunchnumber_MBDNS_scaled{0};

  std::array<long, 120> m_bunchnumber_crossings{0};

  std::array<double, 120> m_bunchnumber_lumi_raw{0.};
  std::array<double, 120> m_bunchnumber_lumi_live{0.};
  std::array<double, 120> m_bunchnumber_lumi_scaled{0.};

  double m_lumi_raw{0.};
  double m_lumi_live{0.};
  double m_lumi_scaled{0.};

  StreamingLumiInfo *m_streaming_lumi_info = nullptr;

};

#endif  // BCOLUMICOUNT_STREAMINGLUMIRECO_H
