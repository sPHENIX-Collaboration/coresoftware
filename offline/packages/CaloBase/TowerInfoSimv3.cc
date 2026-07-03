#include "TowerInfoSimv3.h"

#include "TowerInfo.h"

#include <phool/phool.h>

#include <TSystem.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>

void TowerInfoSimv3::Reset()
{
  TowerInfoSimv1::Reset();
  std::ranges::fill(_waveform,0);
}

void TowerInfoSimv3::set_nsample(int nsample)
{
  if (nsample > 0)
  {
    _waveform.resize(nsample, 0);
    return;
  }
  std::cout << PHWHERE << " invalid number of samples: " << nsample << std::endl;
  gSystem->Exit(1);
  exit(1);
}

int16_t TowerInfoSimv3::get_waveform_value(int index) const
{
  if (index >= 0 && index < get_nsample())
  {
    return _waveform[index];
  }
  return 0;
}

void TowerInfoSimv3::set_waveform_value(int index, int16_t value)
{
  if (index >= 0 && index < get_nsample())
  {
    _waveform[index] = value;
  }
  return;
}

void TowerInfoSimv3::copy_tower(TowerInfo* tower)
{
  TowerInfoSimv1::copy_tower(tower);
  const int nsamples = tower->get_nsample();
  if (nsamples <= 0)
  {
    _waveform.clear();
    return;
  }
  set_nsample(nsamples);
  for (int i = 0; i < nsamples; ++i)
  {
    _waveform[i] = tower->get_waveform_value(i);
  }
  return;
}
