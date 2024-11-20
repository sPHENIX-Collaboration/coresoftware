///////////////////////////////////////////////////////////////
//
// SpinDBContent class
// Author      : D. Loomis, D. Neff (from Y. Fukao PHENIX class)
// Description : Content of spin database
// Created     : 2024-05-12
//
// ERROR_VALUE       : Error value
// NCROSS            : Number of crossing (120)
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

//#include <stdio.h>
//#include <iostream>

class SpinDBContent
{
 public:
  SpinDBContent() { Initialize(); }
  virtual ~SpinDBContent() { ; }
  void Initialize();
  static int GetNCrossing() { return (NCROSS); }
  static int GetErrorValue() { return (ERROR_VALUE); }
  int CheckBunchNumber(int bunch);
  void Print() const;

  int GetRunNumber() { return (runnum); }
  int GetQALevel() { return (qa_level); }
  int GetFillNumber() { return (fillnum); }
  int GetBadRunFlag() { return (badrun); }
  int GetCrossingShift() { return (cross_shift); }

  int GetPolarizationBlue(int bunch, float &value, float &error);
  int GetPolarizationBlue(int bunch, float &value, float &error, float &syserr);
  int GetPolarizationBlue(int bunch, double &value, double &error);
  int GetPolarizationBlue(int bunch, double &value, double &error, double &syserr);
  int GetPolarizationYellow(int bunch, float &value, float &error);
  int GetPolarizationYellow(int bunch, float &value, float &error, float &syserr);
  int GetPolarizationYellow(int bunch, double &value, double &error);
  int GetPolarizationYellow(int bunch, double &value, double &error, double &syserr);
  int GetSpinPatternBlue(int bunch);
  int GetSpinPatternYellow(int bunch);
  long long GetScalerMbdVertexCut(int bunch);
  long long GetScalerMbdNoCut(int bunch);
  long long GetScalerZdcNoCut(int bunch);
  long long GetScaler(int channel, int bunch);
  int GetBadBunchFlag(int bunch);

  void GetAsymBlueForward(float &value, float &error);
  void GetAsymBlueBackward(float &value, float &error);
  void GetAsymYellowForward(float &value, float &error);
  void GetAsymYellowBackward(float &value, float &error);
  void GetPhaseBlueForward(float &value, float &error);
  void GetPhaseBlueBackward(float &value, float &error);
  void GetPhaseYellowForward(float &value, float &error);
  void GetPhaseYellowBackward(float &value, float &error);

  float GetCrossAngle() { return cross_angle; }
  float GetCrossAngleStd() { return cross_angle_std; }
  float GetCrossAngleMin() { return cross_angle_min; }
  float GetCrossAngleMax() { return cross_angle_max; }

  void SetRunNumber(int run)
  {
    runnum = run;
    return;
  }
  void SetQALevel(int qa)
  {
    qa_level = qa;
    return;
  }
  void SetFillNumber(int fill)
  {
    fillnum = fill;
    return;
  }
  void SetBadRunFlag(int flag)
  {
    badrun = flag;
    return;
  }
  void SetCrossingShift(int shift)
  {
    cross_shift = shift;
    return;
  }
  int SetPolarizationBlue(int bunch, float value, float error);
  int SetPolarizationYellow(int bunch, float value, float error);
  int SetPolarizationBlue(int bunch, float value, float error, float syserr);
  int SetPolarizationYellow(int bunch, float value, float error, float syserr);
  int SetSpinPatternBlue(int bunch, int value);
  int SetSpinPatternYellow(int bunch, int value);
  int SetScalerMbdVertexCut(int bunch, long long value);
  int SetScalerMbdNoCut(int bunch, long long value);
  int SetScalerZdcNoCut(int bunch, long long value);
  int SetScaler(int channel, int bunch, long long value);
  int SetBadBunchFlag(int bunch, int value);

  void SetAsymBlueForward(float value, float error);
  void SetAsymBlueBackward(float value, float error);
  void SetAsymYellowForward(float value, float error);
  void SetAsymYellowBackward(float value, float error);
  void SetPhaseBlueForward(float value, float error);
  void SetPhaseBlueBackward(float value, float error);
  void SetPhaseYellowForward(float value, float error);
  void SetPhaseYellowBackward(float value, float error);

  void SetCrossAngle(float value) { cross_angle = value; }
  void SetCrossAngleStd(float value) { cross_angle_std = value; }
  void SetCrossAngleMin(float value) { cross_angle_min = value; }
  void SetCrossAngleMax(float value) { cross_angle_max = value; }


  



 private:
  static const int NCROSS;
  static const int ERROR_VALUE;

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
};

#endif /* USPIN_SPINDBCONTENT_H */
