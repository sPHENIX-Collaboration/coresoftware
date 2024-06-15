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

  tc_x_blue = (float) ERROR_VALUE;
  tc_x_blueerr = (float) ERROR_VALUE;
  tc_y_blue = (float) ERROR_VALUE;
  tc_y_blueerr = (float) ERROR_VALUE;
  tc_x_yellow = (float) ERROR_VALUE;
  tc_x_yellowerr = (float) ERROR_VALUE;
  tc_y_yellow = (float) ERROR_VALUE;
  tc_y_yellowerr = (float) ERROR_VALUE;

  return;
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

  std::cout << (boost::format("Transvers comp. blue X   = %f +- %f\n") % tc_x_blue % tc_x_blueerr).str();
  std::cout << (boost::format("Transvers comp. blue Y   = %f +- %f\n") % tc_y_blue % tc_y_blueerr).str();
  std::cout << (boost::format("Transvers comp. yellow X = %f +- %f\n") % tc_x_yellow % tc_x_yellowerr).str();
  std::cout << (boost::format("Transvers comp. yellow Y = %f +- %f\n") % tc_y_yellow % tc_y_yellowerr).str();

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

void SpinDBContent::SetTransCompBlueX(float value, float error)
{
  tc_x_blue = value;
  tc_x_blueerr = error;
  return;
}

//////////////////////////////////////////////////////

void SpinDBContent::SetTransCompBlueY(float value, float error)
{
  tc_y_blue = value;
  tc_y_blueerr = error;
  return;
}

////////////////////////////////////////////////////////////////

void SpinDBContent::SetTransCompYellowX(float value, float error)
{
  tc_x_yellow = value;
  tc_x_yellowerr = error;
  return;
}

/////////////////////////////////////////////////////////////////

void SpinDBContent::SetTransCompYellowY(float value, float error)
{
  tc_y_yellow = value;
  tc_y_yellowerr = error;
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

void SpinDBContent::GetTransCompBlueX(float &value, float &error)
{
  value = tc_x_blue;
  error = tc_x_blueerr;
  return;
}

///////////////////////////////////////////////////////////////

void SpinDBContent::GetTransCompBlueX(double &value, double &error)
{
  value = (double) tc_x_blue;
  error = (double) tc_x_blueerr;
  return;
}

///////////////////////////////////////////////////////////////

void SpinDBContent::GetTransCompBlueY(float &value, float &error)
{
  value = tc_y_blue;
  error = tc_y_blueerr;
  return;
}

///////////////////////////////////////////////////////////////

void SpinDBContent::GetTransCompBlueY(double &value, double &error)
{
  value = (double) tc_y_blue;
  error = (double) tc_y_blueerr;
  return;
}

///////////////////////////////////////////////////////////////

void SpinDBContent::GetTransCompYellowX(float &value, float &error)
{
  value = tc_x_yellow;
  error = tc_x_yellowerr;
  return;
}

//////////////////////////////////////////////////////////////

void SpinDBContent::GetTransCompYellowX(double &value, double &error)
{
  value = (double) tc_x_yellow;
  error = (double) tc_x_yellowerr;
  return;
}

//////////////////////////////////////////////////////////////

void SpinDBContent::GetTransCompYellowY(float &value, float &error)
{
  value = tc_y_yellow;
  error = tc_y_yellowerr;
  return;
}

//////////////////////////////////////////////////////////////

void SpinDBContent::GetTransCompYellowY(double &value, double &error)
{
  value = (double) tc_y_yellow;
  error = (double) tc_y_yellowerr;
  return;
}

//////////////////////////////////////////////////////////////
