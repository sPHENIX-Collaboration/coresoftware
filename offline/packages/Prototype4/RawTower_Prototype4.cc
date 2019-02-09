#include "RawTower_Prototype4.h"
#include <calobase/RawTowerDefs.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>

#include "PROTOTYPE4_FEM.h"

using namespace std;

    RawTower_Prototype4::RawTower_Prototype4()
  : towerid(~0)
  ,  // initialize all bits on
  energy(0)
  , time(NAN)
  , HBD_channel(-1)
{
  for (int i = 0; i < NSAMPLES; ++i) signal_samples[i] = -9999;
}

RawTower_Prototype4::RawTower_Prototype4(const RawTower& tower)
{
  towerid = (tower.get_id());
  energy = (tower.get_energy());
  time = (tower.get_time());
  HBD_channel = -1;
  for (int i = 0; i < NSAMPLES; ++i) signal_samples[i] = -9999;
}

RawTower_Prototype4::RawTower_Prototype4(RawTowerDefs::keytype id)
  : towerid(id)
  , energy(0)
  , time(NAN)
  , HBD_channel(-1)
{
  for (int i = 0; i < NSAMPLES; ++i) signal_samples[i] = -9999;
}

RawTower_Prototype4::RawTower_Prototype4(const unsigned int icol, const unsigned int irow)
  : towerid(0)
  , energy(0)
  , time(NAN)
  , HBD_channel(-1)
{
  towerid = RawTowerDefs::encode_towerid(RawTowerDefs::NONE, icol, irow);
  for (int i = 0; i < NSAMPLES; ++i) signal_samples[i] = -9999;
}

RawTower_Prototype4::RawTower_Prototype4(const RawTowerDefs::CalorimeterId caloid,
                                         const unsigned int ieta, const unsigned int iphi)
  : towerid(0)
  , energy(0)
  , time(NAN)
  , HBD_channel(-1)
{
  towerid = RawTowerDefs::encode_towerid(caloid, ieta, iphi);
  for (int i = 0; i < NSAMPLES; ++i) signal_samples[i] = -9999;
}

RawTower_Prototype4::~RawTower_Prototype4()
{
}

void RawTower_Prototype4::Reset()
{
  energy = 0;
  time = NAN;
}

int RawTower_Prototype4::isValid() const
{
  return get_energy() != 0;
}

void RawTower_Prototype4::identify(std::ostream& os) const
{
  os << "RawTower_Prototype4: etabin: " << get_bineta() << ", phibin: " << get_binphi()
     << " energy=" << get_energy() << std::endl;
}

void RawTower_Prototype4::set_signal_samples(int i, RawTower_Prototype4::signal_type sig)
{
  assert(i >= 0);
  assert(i < NSAMPLES);
  signal_samples[i] = sig;
}

RawTower_Prototype4::signal_type
RawTower_Prototype4::get_signal_samples(int i) const
{
  assert(i >= 0);
  assert(i < NSAMPLES);
  return signal_samples[i];
}

double
RawTower_Prototype4::get_energy_power_law_exp(int verbosity)
{
  double peak = NAN;
  double peak_sample = NAN;
  double pedstal = NAN;

  vector<double> vec_signal_samples;
  for (int i = 0; i < NSAMPLES; i++)
  {
    vec_signal_samples.push_back(signal_samples[i]);
  }

  PROTOTYPE4_FEM::
      SampleFit_PowerLawExp(vec_signal_samples, peak, peak_sample, pedstal, verbosity);

  return peak;
}

double
RawTower_Prototype4::get_energy_peak_sample(int verbosity)
{
  double peak = NAN;
  double peak_sample = NAN;
  double pedstal = NAN;

  vector<double> vec_signal_samples;
  for (int i = 0; i < NSAMPLES; i++)
  {
    vec_signal_samples.push_back(signal_samples[i]);
  }

  PROTOTYPE4_FEM::
      SampleFit_PeakSample(vec_signal_samples, peak, peak_sample, pedstal, verbosity);

  return peak;
}

double
RawTower_Prototype4::get_energy_power_law_double_exp(int verbosity)
{
  double peak = NAN;
  double peak_sample = NAN;
  double pedstal = NAN;

  vector<double> vec_signal_samples;
  map<int, double> parameters_io;

  for (int i = 0; i < NSAMPLES; i++)
  {
    vec_signal_samples.push_back(signal_samples[i]);
  }

  PROTOTYPE4_FEM::
      SampleFit_PowerLawDoubleExp(vec_signal_samples, peak, peak_sample, pedstal, parameters_io, verbosity);

  return peak;
}
