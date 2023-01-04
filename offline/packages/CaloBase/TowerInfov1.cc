#include "TowerInfov1.h"

TowerInfov1::TowerInfov1(const TowerInfov1 &ti)
   
= default;

void TowerInfov1::Reset()
{
  _time = 0;
  _amplitude = NAN;
}
