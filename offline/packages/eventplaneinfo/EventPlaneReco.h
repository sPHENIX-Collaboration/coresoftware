// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTPLANERECO_H
#define EVENTPLANERECO_H

//===========================================================
/// \author Ejiro Umaka
//===========================================================

#include <fun4all/SubsysReco.h>

#include <string>  // for string
#include <vector>  // for vector

class PHCompositeNode;

class EventPlaneReco : public SubsysReco
{
 public:
  EventPlaneReco(const std::string &name = "EventPlaneReco");
  ~EventPlaneReco() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void ResetMe();
  void set_sepd_epreco(bool sepdEpReco)
  {
    _sepdEpReco = sepdEpReco;
  }
  void set_mbd_epreco(bool mbdEpReco)
  {
    _mbdEpReco = mbdEpReco;
  }
  void set_sEPD_Mip_cut(const float &e)
  {
    _e = e;
  }

 private:
  int CreateNodes(PHCompositeNode *topNode);
  unsigned int m_MaxOrder = 0;

  std::vector<std::vector<double>> south_q;
  std::vector<std::vector<double>> north_q;

  std::vector<std::pair<double, double>> south_Qvec;
  std::vector<std::pair<double, double>> north_Qvec;
  bool _mbdEpReco = false;
  bool _sepdEpReco = false;
  float _e = 6.0;
};

#endif  // EVENTPLANERECO_H
