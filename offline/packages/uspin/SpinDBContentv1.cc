#include "SpinDBContentv1.h"

#include <phool/PHObject.h>

#include <format>
#include <iostream>

//////////////////////////////////////////////////////////

void SpinDBContentv1::InitializeV1()
{
  runnum = GetErrorValue();
  qa_level = GetErrorValue();
  fillnum = GetErrorValue();
  badrun = GetErrorValue();
  cross_shift = GetErrorValue();

  for (int icross = 0; icross < GetNCrossing(); icross++)
  {
    bpol[icross] = (float) GetErrorValue();
    bpolerr[icross] = (float) GetErrorValue();
    bpolsys[icross] = (float) GetErrorValue();
    ypol[icross] = (float) GetErrorValue();
    ypolerr[icross] = (float) GetErrorValue();
    ypolsys[icross] = (float) GetErrorValue();
    bpat[icross] = GetErrorValue();
    ypat[icross] = GetErrorValue();
    scaler_mbd_vtxcut[icross] = (long long) GetErrorValue();
    scaler_mbd_nocut[icross] = (long long) GetErrorValue();
    scaler_zdc_nocut[icross] = (long long) GetErrorValue();
    bad_bunch[icross] = GetErrorValue();
  }

  cross_angle = (float) GetErrorValue();
  cross_angle_std = (float) GetErrorValue();
  cross_angle_min = (float) GetErrorValue();
  cross_angle_max = (float) GetErrorValue();

  asym_bf = (float) GetErrorValue();
  asym_bb = (float) GetErrorValue();
  asym_yf = (float) GetErrorValue();
  asym_yb = (float) GetErrorValue();
  asymerr_bf = (float) GetErrorValue();
  asymerr_bb = (float) GetErrorValue();
  asymerr_yf = (float) GetErrorValue();
  asymerr_yb = (float) GetErrorValue();
  phase_bf = (float) GetErrorValue();
  phase_bb = (float) GetErrorValue();
  phase_yf = (float) GetErrorValue();
  phase_yb = (float) GetErrorValue();
  phaseerr_bf = (float) GetErrorValue();
  phaseerr_bb = (float) GetErrorValue();
  phaseerr_yf = (float) GetErrorValue();
  phaseerr_yb = (float) GetErrorValue();
}

/////////////////////////////////////////////////////////////////

int SpinDBContentv1::CheckBunchNumber(int bunch) const
{
  if (bunch < 0 || bunch >= GetNCrossing())
  {
    std::cout << std::format("Error : bunch number ({}) is out of range (0-119).", bunch) << std::endl;
    return (GetErrorValue());
  }

  return (1);
}

///////////////////////////////////////////////////////////////

void SpinDBContentv1::identify(std::ostream &os) const
{
  os << std::format("Run number = {}", runnum) << std::endl;
  os << std::format("QA Level = {}", qa_level) << std::endl;
  os << std::format("Fill number = {}", fillnum) << std::endl;
  os << std::format("Bad run QA = {}", badrun) << std::endl;
  os << std::format("Crossing shift = {}", cross_shift) << std::endl;

  for (int i = 0; i < GetNCrossing(); i++)
  {
    os << std::format("{:3} : {:12} {:12} {:12} : {:3} {:3} : {}", i, scaler_mbd_vtxcut[i], scaler_mbd_nocut[i], scaler_zdc_nocut[i], bpat[i], ypat[i], bad_bunch[i]) << std::endl;

    os << std::format(" : {:6.3f} +- {:6.3f} +- {:6.3f} {:6.3f} +- {:6.3} +- {:6.3f}", bpol[i], bpolerr[i], bpolsys[i], ypol[i], ypolerr[i], ypolsys[i]) << std::endl;
  }

  return;
}

///////////////////////////////////////////////////////////

int SpinDBContentv1::SetPolarizationBlue(int bunch, float value, float error)
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  bpol[bunch] = value;
  bpolerr[bunch] = error;
  return (1);
}

///////////////////////////////////////////////////////////

int SpinDBContentv1::SetPolarizationBlue(int bunch, float value, float error, float syserr)
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  bpol[bunch] = value;
  bpolerr[bunch] = error;
  bpolsys[bunch] = syserr;
  return (1);
}

/////////////////////////////////////////////////////

int SpinDBContentv1::SetPolarizationYellow(int bunch, float value, float error)
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  ypol[bunch] = value;
  ypolerr[bunch] = error;
  return (1);
}
/////////////////////////////////////////////////////

int SpinDBContentv1::SetPolarizationYellow(int bunch, float value, float error, float syserr)
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  ypol[bunch] = value;
  ypolerr[bunch] = error;
  ypolsys[bunch] = syserr;
  return (1);
}
/////////////////////////////////////////////////////

int SpinDBContentv1::SetSpinPatternBlue(int bunch, int value)
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  bpat[bunch] = value;
  return (1);
}

/////////////////////////////////////////////////////

int SpinDBContentv1::SetSpinPatternYellow(int bunch, int value)
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  ypat[bunch] = value;
  return (1);
}

/////////////////////////////////////////////////////

int SpinDBContentv1::SetScalerMbdVertexCut(int bunch, long long value)
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  scaler_mbd_vtxcut[bunch] = value;
  return (1);
}

//////////////////////////////////////////////////////

int SpinDBContentv1::SetScalerMbdNoCut(int bunch, long long value)
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  scaler_mbd_nocut[bunch] = value;
  return (1);
}

//////////////////////////////////////////////////////

int SpinDBContentv1::SetScalerZdcNoCut(int bunch, long long value)
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  scaler_zdc_nocut[bunch] = value;
  return (1);
}

////////////////////////////////////////////////////////

int SpinDBContentv1::SetScaler(int channel, int bunch, long long value)
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
  return GetErrorValue();
}

////////////////////////////////////////////////////////////////

int SpinDBContentv1::SetBadBunchFlag(int bunch, int value)
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  bad_bunch[bunch] = value;
  return (1);
}

//////////////////////////////////////////////////////

void SpinDBContentv1::SetAsymBlueForward(float value, float error)
{
  asym_bf = value;
  asymerr_bf = error;
  return;
}

//////////////////////////////////////////////////////

void SpinDBContentv1::SetAsymBlueBackward(float value, float error)
{
  asym_bb = value;
  asymerr_bb = error;
  return;
}

////////////////////////////////////////////////////////////////

void SpinDBContentv1::SetAsymYellowForward(float value, float error)
{
  asym_yf = value;
  asymerr_yf = error;
  return;
}

/////////////////////////////////////////////////////////////////

void SpinDBContentv1::SetAsymYellowBackward(float value, float error)
{
  asym_yb = value;
  asymerr_yb = error;
  return;
}

//////////////////////////////////////////////////////

void SpinDBContentv1::SetPhaseBlueForward(float value, float error)
{
  phase_bf = value;
  phaseerr_bf = error;
  return;
}

//////////////////////////////////////////////////////

void SpinDBContentv1::SetPhaseBlueBackward(float value, float error)
{
  phase_bb = value;
  phaseerr_bb = error;
  return;
}

////////////////////////////////////////////////////////////////

void SpinDBContentv1::SetPhaseYellowForward(float value, float error)
{
  phase_yf = value;
  phaseerr_yf = error;
  return;
}

/////////////////////////////////////////////////////////////////

void SpinDBContentv1::SetPhaseYellowBackward(float value, float error)
{
  phase_yb = value;
  phaseerr_yb = error;
  return;
}

////////////////////////////////////////////////////////////////

int SpinDBContentv1::GetPolarizationBlue(int bunch, float &value, float &error) const
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  value = bpol[bunch];
  error = bpolerr[bunch];
  return (1);
}

////////////////////////////////////////////////////////////////

int SpinDBContentv1::GetPolarizationBlue(int bunch, float &value, float &error, float &syserr) const
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  value = bpol[bunch];
  error = bpolerr[bunch];
  syserr = bpolsys[bunch];
  return (1);
}

/////////////////////////////////////////////////////////////////

int SpinDBContentv1::GetPolarizationBlue(int bunch, double &value, double &error) const
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  value = (double) bpol[bunch];
  error = (double) bpolerr[bunch];
  return (1);
}

/////////////////////////////////////////////////////////////////

int SpinDBContentv1::GetPolarizationBlue(int bunch, double &value, double &error, double &syserr) const
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  value = (double) bpol[bunch];
  error = (double) bpolerr[bunch];
  syserr = (double) bpolsys[bunch];
  return (1);
}

/////////////////////////////////////////////////////////////////

int SpinDBContentv1::GetPolarizationYellow(int bunch, float &value, float &error) const
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  value = ypol[bunch];
  error = ypolerr[bunch];
  return (1);
}

/////////////////////////////////////////////////////////////////

int SpinDBContentv1::GetPolarizationYellow(int bunch, float &value, float &error, float &syserr) const
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  value = ypol[bunch];
  error = ypolerr[bunch];
  syserr = ypolsys[bunch];
  return (1);
}

/////////////////////////////////////////////////////////////////

int SpinDBContentv1::GetPolarizationYellow(int bunch, double &value, double &error) const
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  value = (double) ypol[bunch];
  error = (double) ypolerr[bunch];
  return (1);
}

/////////////////////////////////////////////////////////////////

int SpinDBContentv1::GetPolarizationYellow(int bunch, double &value, double &error, double &syserr) const
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  value = (double) ypol[bunch];
  error = (double) ypolerr[bunch];
  syserr = (double) ypolsys[bunch];
  return (1);
}

/////////////////////////////////////////////////////////////////

int SpinDBContentv1::GetSpinPatternBlue(int bunch) const
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  return (bpat[bunch]);
}

////////////////////////////////////////////////////////////////

int SpinDBContentv1::GetSpinPatternYellow(int bunch) const
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  return (ypat[bunch]);
}

////////////////////////////////////////////////////////////////

long long SpinDBContentv1::GetScalerMbdVertexCut(int bunch) const
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  return (scaler_mbd_vtxcut[bunch]);
}

////////////////////////////////////////////////////////////////

long long SpinDBContentv1::GetScalerMbdNoCut(int bunch) const
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  return (scaler_mbd_nocut[bunch]);
}

////////////////////////////////////////////////////////////////

long long SpinDBContentv1::GetScalerZdcNoCut(int bunch) const
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  return (scaler_zdc_nocut[bunch]);
}

////////////////////////////////////////////////////////////////

long long SpinDBContentv1::GetScaler(int channel, int bunch) const
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
  return GetErrorValue();
}

////////////////////////////////////////////////////////////////

int SpinDBContentv1::GetBadBunchFlag(int bunch) const
{
  if (CheckBunchNumber(bunch) == GetErrorValue())
  {
    return (GetErrorValue());
  }
  return (bad_bunch[bunch]);
}

///////////////////////////////////////////////////////////////

void SpinDBContentv1::GetAsymBlueForward(float &value, float &error) const
{
  value = asym_bf;
  error = asymerr_bf;
  return;
}

///////////////////////////////////////////////////////////////

void SpinDBContentv1::GetAsymBlueBackward(float &value, float &error) const
{
  value = asym_bb;
  error = asymerr_bb;
  return;
}

///////////////////////////////////////////////////////////////

void SpinDBContentv1::GetAsymYellowForward(float &value, float &error) const
{
  value = asym_yf;
  error = asymerr_yf;
  return;
}

//////////////////////////////////////////////////////////////

void SpinDBContentv1::GetAsymYellowBackward(float &value, float &error) const
{
  value = asym_yb;
  error = asymerr_yb;
  return;
}

//////////////////////////////////////////////////////////////

void SpinDBContentv1::GetPhaseBlueForward(float &value, float &error) const
{
  value = phase_bf;
  error = phaseerr_bf;
  return;
}

///////////////////////////////////////////////////////////////

void SpinDBContentv1::GetPhaseBlueBackward(float &value, float &error) const
{
  value = phase_bb;
  error = phaseerr_bb;
  return;
}

///////////////////////////////////////////////////////////////

void SpinDBContentv1::GetPhaseYellowForward(float &value, float &error) const
{
  value = phase_yf;
  error = phaseerr_yf;
  return;
}

//////////////////////////////////////////////////////////////

void SpinDBContentv1::GetPhaseYellowBackward(float &value, float &error) const
{
  value = phase_yb;
  error = phaseerr_yb;
  return;
}

//////////////////////////////////////////////////////////////
