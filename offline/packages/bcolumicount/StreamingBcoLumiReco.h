#ifndef BCOLUMICOUNT_STREAMINGBCOLUMIRECO_H
#define BCOLUMICOUNT_STREAMINGBCOLUMIRECO_H

#include "StreamingLumiInfo.h"

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllHistoManager.h>

#include <string>
#include <utility>
#include <array>

class TH1;

class StreamingBcoLumiReco : public SubsysReco
{
 public:
  StreamingBcoLumiReco(const std::string &name = "STREAMINGBCOLUMIRECO");
  ~StreamingBcoLumiReco() override = default;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;

  virtual uint64_t get_bco() const { return m_bco; }

  virtual int get_evtno() const { return m_evtno; }

  virtual bool get_usable_bco_tag() const { return m_usable_bco_tag; }

  virtual std::pair<uint64_t, uint64_t> get_bco_streaming_window() const { return m_bco_streaming_window; }

  virtual double get_lumi_raw() const { return m_lumi_raw; }
  virtual double get_lumi_live() const { return m_lumi_live; }
  virtual double get_lumi_scaled() const { return m_lumi_scaled; }

  virtual void set_default_positive_window_length(int val) { m_default_positive_window_length = val; }
  virtual void set_default_negative_window_length(int val) { m_default_negative_window_length = val; }



 private:
  static int CreateNodeTree(PHCompositeNode *topNode);
  const int trigbits = 40;
  Fun4AllHistoManager *hm = nullptr; 
  TH1 *h_bco_diff = nullptr;
  TH1 *h_bco_diff_trigbits[40] = {nullptr};
  TH1 *h_bco_tag = nullptr;

  uint64_t m_bco{0};
  int m_bunches = 120;
  int m_evtno{0};
  bool m_usable_bco_tag = false;
  std::pair<uint64_t, uint64_t> m_bco_streaming_window;
  unsigned int m_default_positive_window_length{340};
  unsigned int m_default_negative_window_length{20};

  double m_xsec_MBDNS = 24.07*1e9; //convert to pb from Vernier scan DOUBLE CHECK VALUE!

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

#endif  // BCOLUMICOUNT_STREAMINGBCOLUMIRECO_H
