// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CENTRALITYCALIB_MBDCENTRALITYTREE_H
#define CENTRALITYCALIB_MBDCENTRALITYTREE_H

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;

// Writes the per-event MBD/ZDC tree the Au+Au centrality calibration (CentralityCalibrator) runs on: one entry per
// event with any of the MBD trigger bits set (scaled GL1 vector), holding the MBD vertex and t0, the per-tube charge
// and time, and the ZDC energy sums. Needs MbdReco, ZdcReco and GlobalVertexReco upstream.
class MbdCentralityTree : public SubsysReco
{
 public:
  explicit MbdCentralityTree(const std::string &name = "MbdCentralityTree", const std::string &filename = "mbd_centrality_tree.root");
  ~MbdCentralityTree() override = default;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  //! scaled GL1 bits of which at least one must be set for an event to be stored (default 10, 12, 13, 14)
  void SetTriggerBits(const std::vector<int> &bits) { m_trigger_bits = bits; }

 private:
  std::string m_filename;
  std::vector<int> m_trigger_bits{10, 12, 13, 14};

  TFile *m_file{nullptr};
  TTree *m_tree{nullptr};
  long m_nevents{0};

  uint64_t m_gl1_scaled{0};
  float m_mbd_vertex{0};
  float m_mbd_time_zero{0};
  float m_mbd_charge[128]{};
  float m_mbd_time[128]{};
  float m_zdc_sum[2]{};
};

#endif
