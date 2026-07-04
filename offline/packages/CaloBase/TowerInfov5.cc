#include "TowerInfov5.h"
#include "TowerInfo.h"

#include <phool/phool.h>

#include <TSystem.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>

void TowerInfov5::Reset()
{
  TowerInfov2::Reset();
  std::ranges::fill(_waveform, 0);
}

void TowerInfov5::set_nsample(int nsample)
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

int16_t TowerInfov5::get_waveform_value(int index) const
{
  if (index >= 0 && index < get_nsample())
  {
    return _waveform[index];
  }
  return 0;
}

void TowerInfov5::set_waveform_value(int index, int16_t value)
{
  if (index >= 0 && index < get_nsample())
  {
    _waveform[index] = value;
  }
  return;
}

void TowerInfov5::copy_tower(TowerInfo* tower)
{
  TowerInfov2::copy_tower(tower);
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

void TowerInfov5::identify(std::ostream& os) const
{
  os << "TowerInfov5" << std::endl;
  for (int i = 0; i < get_nsample(); ++i)
  {
    os << "sample " << i << ": " << get_waveform_value(i) << std::endl;
  }
  return;
}
