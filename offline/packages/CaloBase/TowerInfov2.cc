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
