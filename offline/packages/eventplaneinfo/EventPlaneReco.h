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

class CDBHistos;
class TProfile2D;
class TH1;

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
  void set_isSim(bool isSim) { _isSim = isSim; }
  void set_sEPD_Mip_cut(const float e) { _epd_e = e; }
  void set_sEPD_Charge_cut(const float c) { _epd_charge_min = c; }
  void set_MBD_Min_Qcut(const float f) { _mbd_e = f; }
  void set_MBD_Vetex_cut(const float v) { _mbd_vertex_cut = v; }
  void set_Ep_orders(const unsigned int n) { m_MaxOrder = n; }

private:
  int CreateNodes(PHCompositeNode *topNode);
  unsigned int m_MaxOrder{3};
  std::string FileName;
    
  std::vector<std::vector<double>> south_q;
  std::vector<std::vector<double>> north_q;
  std::vector<std::vector<double>> northsouth_q;
  static const int nRings = 16;
  std::vector<std::vector<std::vector<double>>> ring_q_north;
  std::vector<std::vector<std::vector<double>>> ring_q_south;
  std::vector<std::pair<double, double>> south_Qvec;
  std::vector<std::pair<double, double>> north_Qvec;
  std::vector<std::pair<double, double>> northsouth_Qvec;
  std::vector<std::vector<std::pair<double, double>>> all_ring_Qvecs_north;
  std::vector<std::vector<std::pair<double, double>>> all_ring_Qvecs_south;
    
  //  const int phibins{24};
  TH1* h_phi_weight_south_input{nullptr};
  TH1* h_phi_weight_north_input{nullptr};
    
  // recentering utility
  std::vector<std::vector<double>> south_q_subtract;
  std::vector<std::vector<double>> north_q_subtract;
  std::vector<std::vector<double>> northsouth_q_subtract;

  // shifting utility
  std::vector<double> shift_north;
  std::vector<double> shift_south;
  std::vector<double> shift_northsouth;
  std::vector<double> tmp_south_psi;
  std::vector<double> tmp_north_psi;
  std::vector<double> tmp_northsouth_psi;

  // recentering histograms
  TProfile2D *tprof_mean_cos_north_epd_input[6]{};
  TProfile2D *tprof_mean_sin_north_epd_input[6]{};
  TProfile2D *tprof_mean_cos_south_epd_input[6]{};
  TProfile2D *tprof_mean_sin_south_epd_input[6]{};
  TProfile2D *tprof_mean_cos_northsouth_epd_input[6]{};
  TProfile2D *tprof_mean_sin_northsouth_epd_input[6]{};

  // shifting histograms
  const int _imax{12};
  TProfile2D *tprof_cos_north_epd_shift_input[6][12]{};
  TProfile2D *tprof_sin_north_epd_shift_input[6][12]{};
  TProfile2D *tprof_cos_south_epd_shift_input[6][12]{};
  TProfile2D *tprof_sin_south_epd_shift_input[6][12]{};
  TProfile2D *tprof_cos_northsouth_epd_shift_input[6][12]{};
  TProfile2D *tprof_sin_northsouth_epd_shift_input[6][12]{};

  bool _mbdEpReco{false};
  bool _sepdEpReco{false};
  bool _isSim{false};
  bool _do_ep{false};

  float _nsum{0.0};
  float _ssum{0.0};
  float _mbdvtx{999.0};
  float _epd_charge_min{5.0};
  float _epd_charge_max{10000.0};
  float _epd_e{10.0};
  float _mbd_e{10.0};
  float _mbdQ{0.0};
  double _totalcharge{0.0};
  float _mbd_vertex_cut{60.0};
};

#endif // EVENTPLANERECO_H
