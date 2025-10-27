#ifndef TIMINGCUT_H
#define TIMINGCUT_H

#include <fun4all/SubsysReco.h>
#include <globalvertex/GlobalVertex.h>
#include <phparameter/PHParameters.h>
#include <cmath>
#include <string>
class PHCompositeNode;

class TimingCut : public SubsysReco
{
 public:
  explicit TimingCut(const std::string &jetNodeName, const std::string &name = "TimingCutModule", int debug = 0, bool doAbort = false, GlobalVertex::VTXTYPE vtxtype = GlobalVertex::MBD);

  ~TimingCut() override = default;

  bool Fails_Delta_t(float lead_time, float sub_time)
  {
    return abs(lead_time - sub_time) > _dt_width;
  }

  bool Fails_Lead_t(float lead_time)
  {
    return abs(lead_time + _t_shift) > _t_width;
  }

  bool Fails_Mbd_dt(float lead_time, float mbd_time)
  {
    return abs(lead_time - mbd_time) > _mbd_dt_width;
  }

  void set_t_shift(float new_shift) { _t_shift = new_shift; }
  float get_t_shift() { return _t_shift; }

  void set_t_width(float new_t_width) { _t_width = new_t_width; }
  float get_t_width() { return _t_width; }
  
  void set_dt_width(float new_dt_width) { _dt_width = new_dt_width; }
  float get_dt_width() { return _dt_width; }

  void set_mbd_dt_width(float new_mbd_dt_width) { _mbd_dt_width = new_mbd_dt_width; }
  float get_mbd_dt_width() { return _mbd_dt_width; }

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  int CreateNodeTree(PHCompositeNode *topNode);

  void SetDefaultParams()
  {
    _cutParams.set_int_param("failsLeadtCut", 0);
    _cutParams.set_int_param("failsDeltatCut", 0);
    _cutParams.set_int_param("failsMbdDtCut", 0);
    _cutParams.set_int_param("failsAnyTimeCut",0);
  }
 private:
  bool _doAbort;
  std::string _name;
  int _debug;
  bool _missingInfoWarningPrinted = false;
  std::string _jetNodeName;
  GlobalVertex::VTXTYPE _vtxtype;
  PHParameters _cutParams;
  float _t_width{8.0};
  float _dt_width{3.0};
  float _t_shift{4.0};
  float _mbd_dt_width{3.0};
};

#endif
