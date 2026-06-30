#include "TowerInfoSimv3.h"

#include "TowerInfo.h"

void TowerInfoSimv3::Reset()
{
  TowerInfoSimv1::Reset();
  for (short& i : _waveform)
  {
    i = 0;
  }
}

void TowerInfoSimv3::Clear(Option_t* /*unused*/)
{
  TowerInfoSimv1::Clear();
  for (short& i : _waveform)
  {
    i = 0;
  }
}

void TowerInfoSimv3::set_nsample(int nsample)
{
  if (nsample >= 0)
  {
    _waveform.resize(nsample, 0);
  }
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
    set_nsample(0);
    return;
  }
  set_nsample(nsamples);
  for (int i = 0; i < nsamples; ++i)
  {
    _waveform[i] = tower->get_waveform_value(i);
  }
  return;
}
