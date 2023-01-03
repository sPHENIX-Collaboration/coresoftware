#include "TowerInfov1.h"

TowerInfov1::TowerInfov1(const TowerInfov1 &ti)
  : TowerInfo(ti)
  , _time(ti._time)
  , _amplitude(ti._amplitude)
{
}

void TowerInfov1::Reset()
{
  _time = 0;
  _amplitude = NAN;
}
