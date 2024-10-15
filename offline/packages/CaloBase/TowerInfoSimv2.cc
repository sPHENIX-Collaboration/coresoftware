#include "TowerInfoSimv2.h"

void TowerInfoSimv2::Reset()
{
  TowerInfoSimv1::Reset();
  for (short& i : _waveform)
  {
    i = 0;
  }
}

void TowerInfoSimv2::Clear(Option_t* /*unused*/)
{
  TowerInfoSimv1::Clear();
  for (short& i : _waveform)
  {
    i = 0;
  }
}

int16_t TowerInfoSimv2::get_waveform_value(int index) const
{
  if (index >= 0 && index < nsample)
  {
    return _waveform[index];
  }
  return 0;
}

void TowerInfoSimv2::set_waveform_value(int index, int16_t value)
{
  if (index >= 0 && index < nsample)
  {
    _waveform[index] = value;
  }
  return;
}

void TowerInfoSimv2::copy_tower(TowerInfo* tower)
{
  TowerInfoSimv1::copy_tower(tower);
  for (int i = 0; i < nsample; ++i)
  {
    set_waveform_value(i, tower->get_waveform_value(i));
  }
  return;
}



