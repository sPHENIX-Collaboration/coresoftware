///////////////////////////////////////////////////////////////
//
// SpinDBContent class
// Author      : D. Loomis (from Y. Fukao PHENIX class)
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
// bpat              : Spin pattern at IP12. (NOT at PHENIX)
// ypat              : Spin pattern at IP12. (NOT at PHENIX)
// scaler_mbd_vtxcut : scaler (MBD VTX)
// scaler_mbd_nocut  : scaler (MBD NS)
// scaler_zdc_nocut   : scaler (ZDC NS)
// bad_bunch         : Bad bunch QA, 0:good else:bad
// tc_x_blue         : Transverse component, horizontal direction, blue beam
// tc_x_blueerr      : Error of tc_x_blue
// tc_y_blue         : Transverse component, vertical direction, blue beam
// tc_y_blueerr      : Error of tc_y_blue
// tc_x_yellow       : Transverse component, horizontal direction, yellow beam
// tc_x_yellowerr    : Error of tc_x_yellow
// tc_y_yellow       : Transverse component, vertical direction, yellow beam
// tc_y_yellowerr    : Error of tc_y_yellow

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

  // bunch xing id from SpinDataEventOut::SpinGL1CrossingID
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

  void GetTransCompBlueX(float &value, float &error);
  void GetTransCompBlueX(double &value, double &error);
  void GetTransCompBlueY(float &value, float &error);
  void GetTransCompBlueY(double &value, double &error);
  void GetTransCompYellowX(float &value, float &error);
  void GetTransCompYellowX(double &value, double &error);
  void GetTransCompYellowY(float &value, float &error);
  void GetTransCompYellowY(double &value, double &error);

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
  void SetTransCompBlueX(float value, float error);
  void SetTransCompBlueY(float value, float error);
  void SetTransCompYellowX(float value, float error);
  void SetTransCompYellowY(float value, float error);

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
  float tc_x_blue;
  float tc_x_blueerr;
  float tc_y_blue;
  float tc_y_blueerr;
  float tc_x_yellow;
  float tc_x_yellowerr;
  float tc_y_yellow;
  float tc_y_yellowerr;
};

#endif /* USPIN_SPINDBCONTENT_H */
