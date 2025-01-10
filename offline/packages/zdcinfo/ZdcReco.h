// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ZDCINFO_ZDCRECO_H
#define ZDCINFO_ZDCRECO_H

//===========================================================
/// \author Ejiro Umaka
//===========================================================

#include <cdbobjects/CDBTTree.h>

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_const_cgsm.h>

#include <string>  // for string
#include <vector>  // for vector

class PHCompositeNode;
class Zdcinfo;

class ZdcReco : public SubsysReco
{
 public:
  ZdcReco(const std::string &name = "ZdcReco");
  ~ZdcReco() override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void ResetMe();
  void set_zdc1_cut(const float f) { _zdc1_e = f; }
  void set_zdc2_cut(const float g) { _zdc2_e = g; }
  void set_smd_threshold(const float d) { _smd_e = d; }

 private:
  void CompSmdPos();
  CDBTTree *cdbttree{nullptr};
  Zdcinfo *m_zdcinfo{nullptr};
  std::string m_Detector{"ZDC"};
  std::string m_fieldname;
  std::string m_calibName;
  bool m_overrideCalibName{false};
  bool m_overrideFieldName{false};
  bool smd_north_fired{false};
  bool smd_south_fired{false};
  std::vector<float> vsmdadc;
  std::vector<float> vsmdtime;
  std::vector<float> vzdcadc;
  std::vector<float> vzdctime;
  float smd_adc[32]{0.0};
  float smd_pos[4]{0.0};
  float _radius_south{0.};
  float _radius_north{0.};
  const double _t{17.623};                                // convert to ns
  const double _c{GSL_CONST_CGSM_SPEED_OF_LIGHT / 1.e9};  // speed of light in cm/ns
  float _z_vertex{0.};
  float _zdc1_e{65.0};
  float _zdc2_e{20.0};
  float _smd_e{5.0};
  float _sumS{0.};
  float _sumN{0.};
  float _sumSden{0.};
  float _sumNden{0.};
  float _sumSt{0.};
  float _sumNt{0.};
  int _nhor{0};
  int _nver{0};
  int _shor{0};
  int _sver{0};
};

#endif  // ZDCRECO_H
