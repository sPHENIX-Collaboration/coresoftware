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

#ifndef USPIN_SPINDBCONTENT_H
#define USPIN_SPINDBCONTENT_H

#include <phool/PHObject.h>

#include <limits>

class SpinDBContent : public PHObject
{
 public:
  SpinDBContent() = default;
  virtual ~SpinDBContent() override = default;

  void identify(std::ostream& os = std::cout) const override;

  static constexpr int GetNCrossing() { return 120; }
  static constexpr int GetErrorValue() { return -999; }

  virtual int CheckBunchNumber(int) const { return 0; }

  virtual int GetRunNumber() const = 0;
  virtual int GetQALevel() const = 0;
  virtual int GetFillNumber() const = 0;
  virtual int GetBadRunFlag() const = 0;
  virtual int GetCrossingShift() const = 0;

  virtual int GetPolarizationBlue(int, float&, float&) const { return -1; }
  virtual int GetPolarizationBlue(int, float&, float&, float&) const { return -1; }
  virtual int GetPolarizationBlue(int, double&, double&) const { return -1; }
  virtual int GetPolarizationBlue(int, double&, double&, double&) const { return -1; }
  virtual int GetPolarizationYellow(int, float&, float&) const { return -1; }
  virtual int GetPolarizationYellow(int, float&, float&, float&) const { return -1; }
  virtual int GetPolarizationYellow(int, double&, double&) const { return -1; }
  virtual int GetPolarizationYellow(int, double&, double&, double&) const { return -1; }

  virtual int GetSpinPatternBlue(int) const { return -1; }
  virtual int GetSpinPatternYellow(int) const { return -1; }
  virtual long long GetScalerMbdVertexCut(int) const { return 0; }
  virtual long long GetScalerMbdNoCut(int) const { return 0; }
  virtual long long GetScalerZdcNoCut(int) const { return 0; }
  virtual long long GetScaler(int, int) const { return 0; }
  virtual int GetBadBunchFlag(int) const { return 0; }

  virtual void GetAsymBlueForward(float& v, float& e) const
  {
    v = std::numeric_limits<float>::quiet_NaN();
    e = std::numeric_limits<float>::quiet_NaN();
  }
  virtual void GetAsymBlueBackward(float& v, float& e) const
  {
    v = std::numeric_limits<float>::quiet_NaN();
    e = std::numeric_limits<float>::quiet_NaN();
  }
  virtual void GetAsymYellowForward(float& v, float& e) const
  {
    v = std::numeric_limits<float>::quiet_NaN();
    e = std::numeric_limits<float>::quiet_NaN();
  }
  virtual void GetAsymYellowBackward(float& v, float& e) const
  {
    v = std::numeric_limits<float>::quiet_NaN();
    e = std::numeric_limits<float>::quiet_NaN();
  }
  virtual void GetPhaseBlueForward(float& v, float& e) const
  {
    v = std::numeric_limits<float>::quiet_NaN();
    e = std::numeric_limits<float>::quiet_NaN();
  }
  virtual void GetPhaseBlueBackward(float& v, float& e) const
  {
    v = std::numeric_limits<float>::quiet_NaN();
    e = std::numeric_limits<float>::quiet_NaN();
  }
  virtual void GetPhaseYellowForward(float& v, float& e) const
  {
    v = std::numeric_limits<float>::quiet_NaN();
    e = std::numeric_limits<float>::quiet_NaN();
  }
  virtual void GetPhaseYellowBackward(float& v, float& e) const
  {
    v = std::numeric_limits<float>::quiet_NaN();
    e = std::numeric_limits<float>::quiet_NaN();
  }

  virtual float GetCrossAngle() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual float GetCrossAngleStd() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual float GetCrossAngleMin() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual float GetCrossAngleMax() const { return std::numeric_limits<float>::quiet_NaN(); }

  virtual void SetRunNumber(int) {}
  virtual void SetQALevel(int) {}
  virtual void SetFillNumber(int) {}
  virtual void SetBadRunFlag(int) {}
  virtual void SetCrossingShift(int) {}

  virtual int SetPolarizationBlue(int, float, float) { return -1; }
  virtual int SetPolarizationYellow(int, float, float) { return -1; }
  virtual int SetPolarizationBlue(int, float, float, float) { return -1; }
  virtual int SetPolarizationYellow(int, float, float, float) { return -1; }

  virtual int SetSpinPatternBlue(int, int) { return -1; }
  virtual int SetSpinPatternYellow(int, int) { return -1; }

  virtual int SetScalerMbdVertexCut(int, long long) { return -1; }
  virtual int SetScalerMbdNoCut(int, long long) { return -1; }
  virtual int SetScalerZdcNoCut(int, long long) { return -1; }
  virtual int SetScaler(int, int, long long) { return -1; }
  virtual int SetBadBunchFlag(int, int) { return -1; }

  virtual void SetAsymBlueForward(float, float) {}
  virtual void SetAsymBlueBackward(float, float) {}
  virtual void SetAsymYellowForward(float, float) {}
  virtual void SetAsymYellowBackward(float, float) {}
  virtual void SetPhaseBlueForward(float, float) {}
  virtual void SetPhaseBlueBackward(float, float) {}
  virtual void SetPhaseYellowForward(float, float) {}
  virtual void SetPhaseYellowBackward(float, float) {}

  virtual void SetCrossAngle(float) {}
  virtual void SetCrossAngleStd(float) {}
  virtual void SetCrossAngleMin(float) {}
  virtual void SetCrossAngleMax(float) {}

 private:
  ClassDefOverride(SpinDBContent, 1);
};

#endif /* USPIN_SPINDBCONTENT_H */
