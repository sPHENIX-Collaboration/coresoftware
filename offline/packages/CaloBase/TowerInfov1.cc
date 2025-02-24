#include "TowerInfov1.h"

#include <limits>

TowerInfov1::TowerInfov1(TowerInfo& tower)
  : _time(tower.get_time())
  , _energy(tower.get_energy())
{
}

void TowerInfov1::Reset()
{
  _time = 0;
  _energy = std::numeric_limits<float>::signaling_NaN();
}

void TowerInfov1::Clear(Option_t* /*unused*/)
{
  _time = 0;
  _energy = 0;
}

void TowerInfov1::copy_tower(TowerInfo* tower)
{
  set_time(tower->get_time());
  set_energy(tower->get_energy());
}
