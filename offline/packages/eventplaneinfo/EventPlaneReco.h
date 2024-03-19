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
class CDBHistos;
class TProfile;
class TH1D;

class EventPlaneReco : public SubsysReco
{
 public:
  EventPlaneReco(const std::string &name = "EventPlaneReco");
  ~EventPlaneReco() override = default;
  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End (PHCompositeNode * /*topNode*/) override;

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
    _epd_e = e;
  }
  void set_MBD_Min_Qcut(const float &f)
  {
    _mbd_e = f;
  }
  void set_z_vertex_cut(const float &v)
  {
    _vertex_cut = v;
  }
  void set_Ep_orders(const unsigned int &n)
  {
    m_MaxOrder = n;
  }
  void set_run_number(const unsigned int &r)
  {
    m_runNo = r;
  }
 
  private:
    
  int CreateNodes(PHCompositeNode *topNode);
   
  unsigned int m_MaxOrder = 3;

  unsigned  int m_runNo = 21813;
 
  std::string OutFileName;
  
  CDBHistos *cdbhistosOut = nullptr;

  std::vector<std::vector<double>> south_q;
  std::vector<std::vector<double>> north_q;
  std::vector<std::pair<double, double>> south_Qvec;
  std::vector<std::pair<double, double>> north_Qvec;
    
  //recentering utility
  std::vector<std::vector<double>> south_q_subtract;
  std::vector<std::vector<double>> north_q_subtract;
    
  //shifting utility
  std::vector<double> tmp_south_psi;
  std::vector<double> tmp_north_psi;
  std::vector<double> shift_north;
  std::vector<double> shift_south;
    

  bool _mbdEpReco = false;
  bool _sepdEpReco = false;
    
  float _epd_e  = 6.0;
  float _mbd_e = 10.0;
  float _mbd_vertex = 999.0;
  float _vertex_cut = 100.0;

  float mbd_e_south;
  float mbd_e_north;
  float mbdQ;
 
  TH1D * hvertex = {};

  //recentering histograms
  TProfile * tprof_mean_cos_north_mbd[6] = {};
  TProfile * tprof_mean_sin_north_mbd[6] = {};
  TProfile * tprof_mean_cos_south_mbd[6] = {};
  TProfile * tprof_mean_sin_south_mbd[6] = {};
    
  TProfile * tprof_mean_cos_north_mbd_input[6] = {};
  TProfile * tprof_mean_sin_north_mbd_input[6] = {};
  TProfile * tprof_mean_cos_south_mbd_input[6] = {};
  TProfile * tprof_mean_sin_south_mbd_input[6] = {};
    
  //shifting histograms
  const int _imax = 6;
  TProfile * tprof_cos_north_mbd_shift[6][6] = {};
  TProfile * tprof_sin_north_mbd_shift[6][6] = {};
  TProfile * tprof_cos_south_mbd_shift[6][6] = {};
  TProfile * tprof_sin_south_mbd_shift[6][6] = {};
    
  TProfile * tprof_cos_north_mbd_shift_input[6][6] = {};
  TProfile * tprof_sin_north_mbd_shift_input[6][6] = {};
  TProfile * tprof_cos_south_mbd_shift_input[6][6] = {};
  TProfile * tprof_sin_south_mbd_shift_input[6][6] = {};
    
};

#endif  // EVENTPLANERECO_H
