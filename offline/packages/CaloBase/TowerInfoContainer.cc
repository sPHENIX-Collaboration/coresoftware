#include "TowerInfoContainer.h"

TowerInfoContainer::Map DummyTowerMap;


TowerInfoContainer::ConstIter TowerInfoContainer::begin() const
{
  return DummyTowerMap.end();
}

TowerInfoContainer::ConstIter TowerInfoContainer::find(int /*key*/) const
{
  return DummyTowerMap.end();
}

TowerInfoContainer::ConstIter TowerInfoContainer::end() const
{
  return DummyTowerMap.end();
}

TowerInfoContainer::Iter TowerInfoContainer::begin()
{
  return DummyTowerMap.end();
}

TowerInfoContainer::Iter TowerInfoContainer::find(int /*key*/)
{
  return DummyTowerMap.end();
}

TowerInfoContainer::Iter TowerInfoContainer::end()
{
  return DummyTowerMap.end();
}

TowerInfoContainer::Range TowerInfoContainer::getTowers()
{
  return std::make_pair(DummyTowerMap.begin(), DummyTowerMap.begin());
}
