#include "SpinDBContent.h"

#include <boost/format.hpp>

#include <iostream>

const int SpinDBContent::NCROSS = 120;
const int SpinDBContent::ERROR_VALUE = -999;

//////////////////////////////////////////////////////////

void SpinDBContent::Initialize()
{
  runnum = ERROR_VALUE;
  qa_level = ERROR_VALUE;
  fillnum = ERROR_VALUE;
  badrun = ERROR_VALUE;
  cross_shift = ERROR_VALUE;

  for (int icross = 0; icross < NCROSS; icross++)
  {
    bpol[icross] = (float) ERROR_VALUE;
    bpolerr[icross] = (float) ERROR_VALUE;
    bpolsys[icross] = (float) ERROR_VALUE;
    ypol[icross] = (float) ERROR_VALUE;
    ypolerr[icross] = (float) ERROR_VALUE;
    ypolsys[icross] = (float) ERROR_VALUE;
    bpat[icross] = ERROR_VALUE;
    ypat[icross] = ERROR_VALUE;
    scaler_mbd_vtxcut[icross] = (long long) ERROR_VALUE;
    scaler_mbd_nocut[icross] = (long long) ERROR_VALUE;
    scaler_zdc_nocut[icross] = (long long) ERROR_VALUE;
    bad_bunch[icross] = ERROR_VALUE;
  }

  cross_angle = (float) ERROR_VALUE;
  cross_angle_std = (float) ERROR_VALUE;
  cross_angle_min = (float) ERROR_VALUE;
  cross_angle_max = (float) ERROR_VALUE;

  asym_bf = (float) ERROR_VALUE;
  asym_bb = (float) ERROR_VALUE;
  asym_yf = (float) ERROR_VALUE;
  asym_yb = (float) ERROR_VALUE;
  asymerr_bf = (float) ERROR_VALUE;
  asymerr_bb = (float) ERROR_VALUE;
  asymerr_yf = (float) ERROR_VALUE;
  asymerr_yb = (float) ERROR_VALUE;
  phase_bf = (float) ERROR_VALUE;
  phase_bb = (float) ERROR_VALUE;
  phase_yf = (float) ERROR_VALUE;
  phase_yb = (float) ERROR_VALUE;
  phaseerr_bf = (float) ERROR_VALUE;
  phaseerr_bb = (float) ERROR_VALUE;
  phaseerr_yf = (float) ERROR_VALUE;
  phaseerr_yb = (float) ERROR_VALUE;



}

/////////////////////////////////////////////////////////////////

int SpinDBContent::CheckBunchNumber(int bunch)
{
  if (bunch < 0 || bunch >= NCROSS)
  {
    std::cout << (boost::format("Error : bunch number (%d) is out of range (0-119).") % bunch).str();
    return (ERROR_VALUE);
  }

  return (1);
}

///////////////////////////////////////////////////////////////

void SpinDBContent::Print() const
{
  std::cout << (boost::format("Run number = %d\n") % runnum).str();
  std::cout << (boost::format("QA Level = %d\n") % qa_level).str();
  std::cout << (boost::format("Fill number = %d\n") % fillnum).str();
  std::cout << (boost::format("Bad run QA = %d\n") % badrun).str();
  std::cout << (boost::format("Crossing shift = %d\n") % cross_shift).str();

  for (int i = 0; i < NCROSS; i++)
  {
    std::cout << (boost::format("%3d : %12lld %12lld %12lld : %3d %3d : %d") % i % scaler_mbd_vtxcut[i] % scaler_mbd_nocut[i] % scaler_zdc_nocut[i] % bpat[i] % ypat[i] % bad_bunch[i]).str();

    std::cout << (boost::format(" : %6.3f +- %6.3f +- %6.3f %6.3f +- %6.3f +- %6.3f\n") % bpol[i] % bpolerr[i] % bpolsys[i] % ypol[i] % ypolerr[i] % ypolsys[i]).str();
  }

  return;
}

///////////////////////////////////////////////////////////

int SpinDBContent::SetPolarizationBlue(int bunch, float value, float error)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  bpol[bunch] = value;
  bpolerr[bunch] = error;
  return (1);
}

///////////////////////////////////////////////////////////

int SpinDBContent::SetPolarizationBlue(int bunch, float value, float error, float syserr)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  bpol[bunch] = value;
  bpolerr[bunch] = error;
  bpolsys[bunch] = syserr;
  return (1);
}

/////////////////////////////////////////////////////

int SpinDBContent::SetPolarizationYellow(int bunch, float value, float error)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  ypol[bunch] = value;
  ypolerr[bunch] = error;
  return (1);
}
/////////////////////////////////////////////////////

int SpinDBContent::SetPolarizationYellow(int bunch, float value, float error, float syserr)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  ypol[bunch] = value;
  ypolerr[bunch] = error;
  ypolsys[bunch] = syserr;
  return (1);
}
/////////////////////////////////////////////////////

int SpinDBContent::SetSpinPatternBlue(int bunch, int value)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  bpat[bunch] = value;
  return (1);
}

/////////////////////////////////////////////////////

int SpinDBContent::SetSpinPatternYellow(int bunch, int value)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  ypat[bunch] = value;
  return (1);
}

/////////////////////////////////////////////////////

int SpinDBContent::SetScalerMbdVertexCut(int bunch, long long value)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  scaler_mbd_vtxcut[bunch] = value;
  return (1);
}

//////////////////////////////////////////////////////

int SpinDBContent::SetScalerMbdNoCut(int bunch, long long value)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  scaler_mbd_nocut[bunch] = value;
  return (1);
}

//////////////////////////////////////////////////////

int SpinDBContent::SetScalerZdcNoCut(int bunch, long long value)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  scaler_zdc_nocut[bunch] = value;
  return (1);
}

////////////////////////////////////////////////////////

int SpinDBContent::SetScaler(int channel, int bunch, long long value)
{
  switch (channel)
  {
  case 0:
    return SetScalerMbdVertexCut(bunch, value);
  case 1:
    return SetScalerMbdNoCut(bunch, value);
  case 2:
    return SetScalerZdcNoCut(bunch, value);
  default:
    break;
  }
  return ERROR_VALUE;
}

////////////////////////////////////////////////////////////////

int SpinDBContent::SetBadBunchFlag(int bunch, int value)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  bad_bunch[bunch] = value;
  return (1);
}

//////////////////////////////////////////////////////

void SpinDBContent::SetAsymBlueForward(float value, float error)
{
  asym_bf = value;
  asymerr_bf = error;
  return;
}

//////////////////////////////////////////////////////

void SpinDBContent::SetAsymBlueBackward(float value, float error)
{
  asym_bb = value;
  asymerr_bb = error;
  return;
}

////////////////////////////////////////////////////////////////

void SpinDBContent::SetAsymYellowForward(float value, float error)
{
  asym_yf = value;
  asymerr_yf = error;
  return;
}

/////////////////////////////////////////////////////////////////

void SpinDBContent::SetAsymYellowBackward(float value, float error)
{
  asym_yb = value;
  asymerr_yb = error;
  return;
}


//////////////////////////////////////////////////////

void SpinDBContent::SetPhaseBlueForward(float value, float error)
{
  phase_bf = value;
  phaseerr_bf = error;
  return;
}

//////////////////////////////////////////////////////

void SpinDBContent::SetPhaseBlueBackward(float value, float error)
{
  phase_bb = value;
  phaseerr_bb = error;
  return;
}

////////////////////////////////////////////////////////////////

void SpinDBContent::SetPhaseYellowForward(float value, float error)
{
  phase_yf = value;
  phaseerr_yf = error;
  return;
}

/////////////////////////////////////////////////////////////////

void SpinDBContent::SetPhaseYellowBackward(float value, float error)
{
  phase_yb = value;
  phaseerr_yb = error;
  return;
}

////////////////////////////////////////////////////////////////

int SpinDBContent::GetPolarizationBlue(int bunch, float &value, float &error)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  value = bpol[bunch];
  error = bpolerr[bunch];
  return (1);
}

////////////////////////////////////////////////////////////////

int SpinDBContent::GetPolarizationBlue(int bunch, float &value, float &error, float &syserr)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  value = bpol[bunch];
  error = bpolerr[bunch];
  syserr = bpolsys[bunch];
  return (1);
}

/////////////////////////////////////////////////////////////////

int SpinDBContent::GetPolarizationBlue(int bunch, double &value, double &error)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  value = (double) bpol[bunch];
  error = (double) bpolerr[bunch];
  return (1);
}

/////////////////////////////////////////////////////////////////

int SpinDBContent::GetPolarizationBlue(int bunch, double &value, double &error, double &syserr)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  value = (double) bpol[bunch];
  error = (double) bpolerr[bunch];
  syserr = (double) bpolsys[bunch];
  return (1);
}

/////////////////////////////////////////////////////////////////

int SpinDBContent::GetPolarizationYellow(int bunch, float &value, float &error)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  value = ypol[bunch];
  error = ypolerr[bunch];
  return (1);
}

/////////////////////////////////////////////////////////////////

int SpinDBContent::GetPolarizationYellow(int bunch, float &value, float &error, float &syserr)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  value = ypol[bunch];
  error = ypolerr[bunch];
  syserr = ypolsys[bunch];
  return (1);
}

/////////////////////////////////////////////////////////////////

int SpinDBContent::GetPolarizationYellow(int bunch, double &value, double &error)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  value = (double) ypol[bunch];
  error = (double) ypolerr[bunch];
  return (1);
}

/////////////////////////////////////////////////////////////////

int SpinDBContent::GetPolarizationYellow(int bunch, double &value, double &error, double &syserr)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  value = (double) ypol[bunch];
  error = (double) ypolerr[bunch];
  syserr = (double) ypolsys[bunch];
  return (1);
}

/////////////////////////////////////////////////////////////////

int SpinDBContent::GetSpinPatternBlue(int bunch)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  return (bpat[bunch]);
}

////////////////////////////////////////////////////////////////

int SpinDBContent::GetSpinPatternYellow(int bunch)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  return (ypat[bunch]);
}

////////////////////////////////////////////////////////////////

long long SpinDBContent::GetScalerMbdVertexCut(int bunch)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  return (scaler_mbd_vtxcut[bunch]);
}

////////////////////////////////////////////////////////////////

long long SpinDBContent::GetScalerMbdNoCut(int bunch)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  return (scaler_mbd_nocut[bunch]);
}

////////////////////////////////////////////////////////////////

long long SpinDBContent::GetScalerZdcNoCut(int bunch)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  return (scaler_zdc_nocut[bunch]);
}

////////////////////////////////////////////////////////////////

long long SpinDBContent::GetScaler(int channel, int bunch)
{
  switch (channel)
  {
  case 0:
    return GetScalerMbdVertexCut(bunch);
  case 1:
    return GetScalerMbdNoCut(bunch);
  case 2:
    return GetScalerZdcNoCut(bunch);
  default:
    break;
  }
  return ERROR_VALUE;
}

////////////////////////////////////////////////////////////////

int SpinDBContent::GetBadBunchFlag(int bunch)
{
  if (CheckBunchNumber(bunch) == ERROR_VALUE)
  {
    return (ERROR_VALUE);
  }
  return (bad_bunch[bunch]);
}

///////////////////////////////////////////////////////////////

void SpinDBContent::GetAsymBlueForward(float &value, float &error)
{
  value = asym_bf;
  error = asymerr_bf;
  return;
}

///////////////////////////////////////////////////////////////

void SpinDBContent::GetAsymBlueBackward(float &value, float &error)
{
  value = asym_bb;
  error = asymerr_bb;
  return;
}

///////////////////////////////////////////////////////////////

void SpinDBContent::GetAsymYellowForward(float &value, float &error)
{
  value = asym_yf;
  error = asymerr_yf;
  return;
}

//////////////////////////////////////////////////////////////

void SpinDBContent::GetAsymYellowBackward(float &value, float &error)
{
  value = asym_yb;
  error = asymerr_yb;
  return;
}

//////////////////////////////////////////////////////////////

void SpinDBContent::GetPhaseBlueForward(float &value, float &error)
{
  value = phase_bf;
  error = phaseerr_bf;
  return;
}

///////////////////////////////////////////////////////////////

void SpinDBContent::GetPhaseBlueBackward(float &value, float &error)
{
  value = phase_bb;
  error = phaseerr_bb;
  return;
}

///////////////////////////////////////////////////////////////

void SpinDBContent::GetPhaseYellowForward(float &value, float &error)
{
  value = phase_yf;
  error = phaseerr_yf;
  return;
}

//////////////////////////////////////////////////////////////

void SpinDBContent::GetPhaseYellowBackward(float &value, float &error)
{
  value = phase_yb;
  error = phaseerr_yb;
  return;
}

//////////////////////////////////////////////////////////////
