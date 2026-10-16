#ifndef BCOLUMICOUNT_STREAMINGBCORECO_H
#define BCOLUMICOUNT_STREAMINGBCORECO_H

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllHistoManager.h>

#include <string>
#include <utility>
#include <array>

class TH1;

class StreamingBcoReco : public SubsysReco
{
 public:
  StreamingBcoReco(const std::string &name = "STREAMINGBCOLUMIRECO");
  ~StreamingBcoReco() override = default;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  virtual int get_evtno() const { return m_evtno; }

  virtual uint64_t get_bco() const { return m_bco; }

  virtual bool get_usable_bco_tag() const { return m_usable_bco_tag; }

  virtual std::pair<uint64_t, uint64_t> get_bco_streaming_window() const { return m_bco_streaming_window; }

  virtual void set_default_positive_window_length(int val) { m_default_positive_window_length = val; }
  virtual void set_default_negative_window_length(int val) { m_default_negative_window_length = val; }



 private:
  static int CreateNodeTree(PHCompositeNode *topNode);
  Fun4AllHistoManager *hm = nullptr; 
  TH1 *h_bco_diff = nullptr;
  TH1 *h_bco_tag = nullptr;

  uint64_t m_bco{0};
  int m_evtno{0};
  bool m_usable_bco_tag = false;
  std::pair<uint64_t, uint64_t> m_bco_streaming_window;
  // max BCO for streaming window = (983 (min acquisition window) - 254.19 (max drift edge) + 8 (min T0) ) / 30*8 (ADC->BCO) = 392 BCO
  // Values taken from Yuri Mitrankov
  unsigned int m_default_positive_window_length{392};
  unsigned int m_default_negative_window_length{0};
};

#endif  // BCOLUMICOUNT_STREAMINGBCORECO_H
