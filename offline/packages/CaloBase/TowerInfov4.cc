#include "TowerInfov4.h"
#include "TowerInfo.h"

#include <limits>

void TowerInfov4::Reset()
{
  energy = 0;
  time = 0;
  chi2 = 0;
  status = 0;
}

void TowerInfov4::copy_tower(TowerInfo* tower)
{
  set_time(tower->get_time());
  set_energy(tower->get_energy());
  set_chi2(tower->get_chi2());
  set_status(tower->get_status());
  return;
}
