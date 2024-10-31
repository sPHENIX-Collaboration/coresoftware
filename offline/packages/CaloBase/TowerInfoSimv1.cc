#include "TowerInfoSimv1.h"

void TowerInfoSimv1::Reset()
{
  TowerInfov2::Reset();
  _hitedeps.clear();
  _showeredeps.clear();
  return;
}

void TowerInfoSimv1::Clear(Option_t* /*unused*/)
{
  TowerInfov2::Clear();
  _hitedeps.clear();
  _showeredeps.clear();
  return;
}


void TowerInfoSimv1::copy_tower(TowerInfo* tower)
{
  TowerInfov2::copy_tower(tower);
  _hitedeps = tower->get_hitEdepMap();
  _showeredeps = tower->get_showerEdepMap();
  return;
}

TowerInfoSimv1::EdepMap& TowerInfoSimv1::get_hitEdepMap()
{
  return _hitedeps;
}

TowerInfoSimv1::ShowerEdepMap& TowerInfoSimv1::get_showerEdepMap()
{
  return _showeredeps;
}

const TowerInfoSimv1::EdepMap& TowerInfoSimv1::get_hitEdepMap() const
{
  return _hitedeps;
}

const TowerInfoSimv1::ShowerEdepMap& TowerInfoSimv1::get_showerEdepMap() const
{
  return _showeredeps;
}

void TowerInfoSimv1::add_edep(const PHG4HitDefs::keytype g4hitid, const float edep)
{
  _hitedeps[g4hitid] += edep;
  return;
}

void TowerInfoSimv1::add_shower_edep(const int showerid, const float edep)
{
  _showeredeps[showerid] += edep;
  return;
}

