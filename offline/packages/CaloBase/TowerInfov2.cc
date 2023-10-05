#include "TowerInfo.h"
#include "TowerInfov2.h"

void TowerInfov2::Reset()
{
  TowerInfov1::Reset();
  _chi2 = 0;
  _pedestal = 0;
  _status = 0;
}

void TowerInfov2::Clear(Option_t* )
{
  TowerInfov1::Clear();
  _chi2 = 0;
  _pedestal = 0;
  _status = 0;
}

void TowerInfov2::copy_tower(TowerInfo* tower)
{
  TowerInfov1::copy_tower(tower);
  set_time_float(tower->get_time_float());
  set_energy(tower->get_energy());
  set_chi2(tower->get_chi2());
  set_pedestal(tower->get_pedestal());
  set_status(tower->get_status());
  return;
}
