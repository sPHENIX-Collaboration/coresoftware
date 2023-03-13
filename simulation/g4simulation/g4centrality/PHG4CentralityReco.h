#ifndef G4CENTRALITY_PHG4CENTRALITYRECO_H
#define G4CENTRALITY_PHG4CENTRALITYRECO_H

//===========================================================
/// \file PHG4CentralityReco.h
/// \brief Centrality quantity construction & calibration
/// \author Dennis V. Perepelitsa
//===========================================================

#include <fun4all/SubsysReco.h>

#include <phparameter/PHParameters.h>

#include <cmath>
#include <map>
#include <string>

class PHCompositeNode;

class PHG4CentralityReco : public SubsysReco
{
 public:
  PHG4CentralityReco(const std::string &name = "PHG4CentralityReco");
  ~PHG4CentralityReco() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void DoCentralityCalibration(bool do_centrality_calibration)
  {
    _do_centrality_calibration = do_centrality_calibration;
  }

  PHParameters &GetCalibrationParameters()
  {
    return _centrality_calibration_params;
  }

 private:
  void CreateNode(PHCompositeNode *topNode);
  void FillNode(PHCompositeNode *topNode);

  PHParameters _centrality_calibration_params;

  bool _do_centrality_calibration = true;

  std::map<float, int> _cent_cal_bimp;
  std::map<float, int> _cent_cal_epd;
  std::map<float, int> _cent_cal_mbd;

  float _bimp = NAN;
  float _bimp_cent = NAN;

  float _epd_N = NAN;
  float _epd_S = NAN;
  float _epd_NS = NAN;
  float _epd_cent = NAN;

  float _mbd_N = NAN;
  float _mbd_S = NAN;
  float _mbd_NS = NAN;
  float _mbd_cent = NAN;
};

#endif
