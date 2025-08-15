///////////////////////////////////////////////////////////////
//
// SpinDBContent class
// Author      : D. Loomis, D. Neff (from Y. Fukao PHENIX class)
// Description : Content of spin database
// Created     : 2024-05-12
//
// runnum            : Run number
// qa_level          : Level of QA for stored data.
// fillnum           : Fill number
// badrun            : Run QA by spin analysis, 1:bad 0:good
// cross_shift       : Crossing sfhit, corrected=(raw+cross_shift)%NCROSS
// bpol              : Blue beam polarization
// bpolerr           : Blue beam polarization error
// bpolsys           : Blue beam polarization systematic error
// ypol              : Yellow beam polarization
// ypolerr           : Yellow beam polarization error
// ypolsys           : Yellow beam polarization systematic error
// bpat              : Spin pattern at IP12. (NOT at sPHENIX)
// ypat              : Spin pattern at IP12. (NOT at sPHENIX)
// scaler_mbd_vtxcut : scaler (MBD VTX)
// scaler_mbd_nocut  : scaler (MBD NS)
// scaler_zdc_nocut  : scaler (ZDC NS)
// bad_bunch         : Bad bunch QA, 0:good else:bad
// cross_angle       : Average relative crossing angle between blue and yellow beams in mrad for run. Sign is dictated by CAD convention, can be + or -
// cross_angle_std   : Standard deviation of relative crossing angle in mrad
// cross_angle_min   : Minimum value of relative crossing angle in mrad
// cross_angle_max   : Maximum value of relative crossing angle in mrad
// asym_bf           : very forward neutron TSSA magnitude in blue beam
// asym_bb           : very backward neutron TSSA magnitude in blue beam
// asym_yf           : very forward neutron TSSA magnitude in yellow beam
// asym_yb           : very backward neutron TSSA magnitude in yellow beam
// asymerr_bf        : very forward neutron TSSA magnitude uncertainty in blue beam
// asymerr_bb        : very backward neutron TSSA magnitude uncertainty in blue beam
// asymerr_yf        : very forward neutron TSSA magnitude uncertainty in yellow beam
// asymerr_yb        : very backward neutron TSSA magnitude uncertainty in yellow beam
// phase_bf          : very forward neutron TSSA phase offset in blue beam
// phase_bb          : very backward neutron TSSA phase offset in blue beam
// phase_yf          : very forward neutron TSSA phase offset in yellow beam
// phase_yb          : very backward neutron TSSA phase offset in yellow beam
// phaseerr_bf       : very forward neutron TSSA phase offset uncertainty in blue beam
// phaseerr_bb       : very backward neutron TSSA phase offset uncertainty in blue beam
// phaseerr_yf       : very forward neutron TSSA phase offset uncertainty in yellow beam
// phaseerr_yb       : very backward neutron TSSA phase offset uncertainty in yellow beam

////////////////////////////////////////////////////////////////

#ifndef USPIN_SPINDBCONTENTV1_H
#define USPIN_SPINDBCONTENTV1_H

#include "SpinDBContent.h"

// #include <stdio.h>
// #include <iostream>

class SpinDBContentv1 : public SpinDBContent
{
 public:
  SpinDBContentv1() { InitializeV1(); }
  ~SpinDBContentv1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void InitializeV1();

  int CheckBunchNumber(int bunch) const override;

  int GetRunNumber() const override { return runnum; }
  int GetQALevel() const override { return qa_level; }
  int GetFillNumber() const override { return fillnum; }
  int GetBadRunFlag() const override { return badrun; }
  int GetCrossingShift() const override { return cross_shift; }

  int GetPolarizationBlue(int, float&, float&) const override;
  int GetPolarizationBlue(int, float&, float&, float&) const override;
  int GetPolarizationBlue(int, double&, double&) const override;
  int GetPolarizationBlue(int, double&, double&, double&) const override;
  int GetPolarizationYellow(int, float&, float&) const override;
  int GetPolarizationYellow(int, float&, float&, float&) const override;
  int GetPolarizationYellow(int, double&, double&) const override;
  int GetPolarizationYellow(int, double&, double&, double&) const override;

  int GetSpinPatternBlue(int) const override;
  int GetSpinPatternYellow(int) const override;
  long long GetScalerMbdVertexCut(int) const override;
  long long GetScalerMbdNoCut(int) const override;
  long long GetScalerZdcNoCut(int) const override;
  long long GetScaler(int, int) const override;
  int GetBadBunchFlag(int) const override;

  void GetAsymBlueForward(float&, float&) const override;
  void GetAsymBlueBackward(float&, float&) const override;
  void GetAsymYellowForward(float&, float&) const override;
  void GetAsymYellowBackward(float&, float&) const override;
  void GetPhaseBlueForward(float&, float&) const override;
  void GetPhaseBlueBackward(float&, float&) const override;
  void GetPhaseYellowForward(float&, float&) const override;
  void GetPhaseYellowBackward(float&, float&) const override;

  float GetCrossAngle() const override { return cross_angle; }
  float GetCrossAngleStd() const override { return cross_angle_std; }
  float GetCrossAngleMin() const override { return cross_angle_min; }
  float GetCrossAngleMax() const override { return cross_angle_max; }

  void SetRunNumber(int run) override { runnum = run; }
  void SetQALevel(int qa) override { qa_level = qa; }
  void SetFillNumber(int fill) override { fillnum = fill; }
  void SetBadRunFlag(int flag) override { badrun = flag; }
  void SetCrossingShift(int shift) override { cross_shift = shift; }

  int SetPolarizationBlue(int, float, float) override;
  int SetPolarizationYellow(int, float, float) override;
  int SetPolarizationBlue(int, float, float, float) override;
  int SetPolarizationYellow(int, float, float, float) override;

  int SetSpinPatternBlue(int, int) override;
  int SetSpinPatternYellow(int, int) override;
  int SetScalerMbdVertexCut(int, long long) override;
  int SetScalerMbdNoCut(int, long long) override;
  int SetScalerZdcNoCut(int, long long) override;
  int SetScaler(int, int, long long) override;
  int SetBadBunchFlag(int, int) override;

  void SetAsymBlueForward(float value, float error) override;
  void SetAsymBlueBackward(float value, float error) override;
  void SetAsymYellowForward(float value, float error) override;
  void SetAsymYellowBackward(float value, float error) override;
  void SetPhaseBlueForward(float value, float error) override;
  void SetPhaseBlueBackward(float value, float error) override;
  void SetPhaseYellowForward(float value, float error) override;
  void SetPhaseYellowBackward(float value, float error) override;

  void SetCrossAngle(float value) override { cross_angle = value; }
  void SetCrossAngleStd(float value) override { cross_angle_std = value; }
  void SetCrossAngleMin(float value) override { cross_angle_min = value; }
  void SetCrossAngleMax(float value) override { cross_angle_max = value; }

 private:
  int runnum;
  int qa_level;
  int fillnum;
  int badrun;
  int cross_shift;
  float bpol[120];
  float bpolerr[120];
  float bpolsys[120];
  float ypol[120];
  float ypolerr[120];
  float ypolsys[120];
  int bpat[120];
  int ypat[120];
  long long scaler_mbd_vtxcut[120];
  long long scaler_mbd_nocut[120];
  long long scaler_zdc_nocut[120];
  int bad_bunch[120];
  float cross_angle;
  float cross_angle_std;
  float cross_angle_min;
  float cross_angle_max;
  float asym_bf;
  float asym_bb;
  float asym_yf;
  float asym_yb;
  float asymerr_bf;
  float asymerr_bb;
  float asymerr_yf;
  float asymerr_yb;
  float phase_bf;
  float phase_bb;
  float phase_yf;
  float phase_yb;
  float phaseerr_bf;
  float phaseerr_bb;
  float phaseerr_yf;
  float phaseerr_yb;

 private:
  ClassDefOverride(SpinDBContentv1, 1);
};

#endif /* USPIN_SPINDBCONTENTV1_H */
