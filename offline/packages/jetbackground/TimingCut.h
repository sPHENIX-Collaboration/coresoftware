#ifndef TIMINGCUT_H
#define TIMINGCUT_H

#include <jetbase/Jet.h>

#include <fun4all/SubsysReco.h>

#include <phparameter/PHParameters.h>

#include <cmath>
#include <string>

class PHCompositeNode;

class TimingCut : public SubsysReco
{
 public:
  explicit TimingCut(const std::string &jetNodeName, const std::string &name = "TimingCutModule", bool doAbort = false);

  ~TimingCut() override = default;

  float calc_dphi(float maxJetPhi, float subJetPhi)
  {
    float dPhi = std::abs(maxJetPhi - subJetPhi);
    if(dPhi>M_PI) dPhi -= M_PI;
    return dPhi;
  }

  bool Pass_Delta_t(float lead_time, float sub_time, float maxJetPhi, float subJetPhi)
  {
    float dPhi = calc_dphi(maxJetPhi, subJetPhi);
    return (std::abs(lead_time - sub_time) < _dt_width && dPhi > _min_dphi);
  }

  bool Pass_Lead_t(float lead_time)
  {
    return std::abs(lead_time + _t_shift) < _t_width;
  }

  bool Pass_Mbd_dt(float lead_time, float mbd_time)
  {
    return std::abs(lead_time - mbd_time) < _mbd_dt_width;
  }

  void set_t_shift(float new_shift) { _t_shift = new_shift; }
  float get_t_shift() { return _t_shift; }

  void set_t_width(float new_t_width) { _t_width = new_t_width; }
  float get_t_width() { return _t_width; }
  
  void set_dt_width(float new_dt_width) { _dt_width = new_dt_width; }
  float get_dt_width() { return _dt_width; }

  void set_mbd_dt_width(float new_mbd_dt_width) { _mbd_dt_width = new_mbd_dt_width; }
  float get_mbd_dt_width() { return _mbd_dt_width; }

  void set_min_dphi(float new_min_dphi) { _min_dphi = new_min_dphi; }
  float get_min_dphi() { return _min_dphi; }

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  int CreateNodeTree(PHCompositeNode *topNode);

  void SetDefaultParams()
  {
    _cutParams.set_int_param("passLeadtCut", 0);
    _cutParams.set_int_param("passDeltatCut", 0);
    _cutParams.set_int_param("passMbdDtCut", 0);
    _cutParams.set_int_param("failAnyTimeCut",0);
    _cutParams.set_double_param("maxJett",9999);
    _cutParams.set_double_param("subJett",9999);
    _cutParams.set_double_param("mbd_time",9999);
    _cutParams.set_double_param("dPhi",9999);
  }

private:
  bool _doAbort;
  bool _missingInfoWarningPrinted = false;
  std::string _jetNodeName;
  PHParameters _cutParams;
  float _t_width{6.0};
  float _dt_width{3.0};
  float _t_shift{2.0};
  float _mbd_dt_width{3.0};
  float _min_dphi{3*M_PI/4};
};

#endif
