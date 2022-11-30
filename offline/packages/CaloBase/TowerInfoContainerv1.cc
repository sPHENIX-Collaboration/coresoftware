#include "TowerInfoContainerv1.h"
#include "TowerInfov1.h"

#include <phool/PHObject.h>

TowerInfoContainerv1::TowerInfoContainerv1()
{
  _clones = new TClonesArray("TowerInfov1", 50);
  _clones->SetOwner();
  _clones->SetName("TowerInfoContainerv1");
}

TowerInfoContainerv1::~TowerInfoContainerv1()
{
  if (_clones)
  {
    delete _clones;
    _clones = 0;
  }

}

void TowerInfoContainerv1::Reset()
{
  _clones->Clear();
}

void TowerInfoContainerv1::add(TowerInfov1 *ti, int pos)
{
  new ((*_clones)[pos]) TowerInfov1(*ti);
}

TowerInfov1* TowerInfoContainerv1::at(int pos)
{
  return (TowerInfov1*)_clones->At(pos);
}
