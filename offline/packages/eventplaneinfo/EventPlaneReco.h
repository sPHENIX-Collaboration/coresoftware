// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTPLANERECO_H
#define EVENTPLANERECO_H

//===========================================================
/// \author Ejiro Umaka
//===========================================================

#include <cdbobjects/CDBTTree.h>
#include <fun4all/SubsysReco.h>

#include <string> // for string
#include <vector> // for vector

class PHCompositeNode;

class EventPlaneReco : public SubsysReco {
public:
  EventPlaneReco(const std::string &name = "EventPlaneReco");
  ~EventPlaneReco() override = default;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode * /*topNode*/) override;

  void ResetMe();
  void set_sepd_epreco(bool sepdEpReco) { _sepdEpReco = sepdEpReco; }
  void set_mbd_epreco(bool mbdEpReco) { _mbdEpReco = mbdEpReco; }
  void set_sEPD_Mip_cut(const float e) { _epd_e = e; }
  void set_MBD_Min_Qcut(const float f) { _mbd_e = f; }
  void set_Ep_orders(const unsigned int n) { m_MaxOrder = n; }

private:
  int CreateNodes(PHCompositeNode *topNode);

  unsigned int m_MaxOrder{3};

  std::vector<std::vector<double>> south_q;
  std::vector<std::vector<double>> north_q;
  std::vector<std::pair<double, double>> south_Qvec;
  std::vector<std::pair<double, double>> north_Qvec;

  bool _mbdEpReco{false};
  bool _sepdEpReco{false};

  float _epd_e{6.0};
  float _mbd_e{10.0};
  float mbdQ{0.};

  std::string m_sEPDMapName;
  std::string m_sEPDfieldname;
  bool m_overrideSEPDMapName{false};
  bool m_overrideSEPDFieldName{false};
  std::vector<unsigned int> vkey;
  unsigned int key{999};
  CDBTTree *cdbttree{nullptr};
};

#endif // EVENTPLANERECO_H
