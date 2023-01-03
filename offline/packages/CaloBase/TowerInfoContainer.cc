#include "TowerInfoContainer.h"
#include "TowerInfo.h"

TowerInfoContainer::TowerMap DummyTowerMap;

void TowerInfoContainer::Reset()
{
}

void TowerInfoContainer::add(TowerInfo* /*ti*/, int /*pos*/)
{
  return;
}

TowerInfo* TowerInfoContainer::at(int /*pos*/)
{
  return nullptr;
}

TowerInfoContainer::TowerMap TowerInfoContainer::getTowerMap()
{
  return DummyTowerMap;
}

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
