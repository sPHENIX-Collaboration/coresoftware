#include "TowerInfov3.h"

void TowerInfov3::Reset()
{
  TowerInfov2::Reset();  
  for (int i = 0; i < nsample; ++i)
  {
    _waveform[i] = 0;
  }
}

void TowerInfov3::Clear(Option_t *)
{
  TowerInfov2::Clear();  
  for (int i = 0; i < nsample; ++i)
  {
    _waveform[i] = 0;
  }
}

int16_t TowerInfov3::get_waveform_value(int index) const
{
  if (index >= 0 && index < nsample)
  {
    return _waveform[index];
  }
  return 0;
}

void TowerInfov3::set_waveform_value(int index, int16_t value)
{
  if (index >= 0 && index < nsample)
  {
    _waveform[index] = value;
  }
  return;
}

void TowerInfov3::copy_tower(TowerInfov3* tower)
{
  TowerInfov2::copy_tower(tower);
  for (int i = 0; i < nsample; ++i)
  {
    set_waveform_value(i, tower->get_waveform_value(i));
  }
  return;
}
