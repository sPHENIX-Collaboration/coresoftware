#include "TowerInfo.h"
#include "TowerInfov4.h"

void TowerInfov4::Reset()
{
  TowerInfov1::Reset();
  _chi2 = 0;
  _status = 0;
}

void TowerInfov4::Clear(Option_t* )
{
  TowerInfov1::Clear();
  _chi2 = 0;
  _status = 0;
}

void TowerInfov4::copy_tower(TowerInfo* tower)
{
  TowerInfov1::copy_tower(tower);
  set_time_float(tower->get_time_float());
  set_energy(tower->get_energy());
  set_chi2(tower->get_chi2());
  set_status(tower->get_status());
  return;
}
