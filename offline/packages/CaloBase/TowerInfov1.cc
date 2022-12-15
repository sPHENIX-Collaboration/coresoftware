#include "TowerInfov1.h"

TowerInfov1::TowerInfov1()
  : _time(0)
  , _amplitude(0)
{
}

TowerInfov1::TowerInfov1(const TowerInfov1 &ti)
  : TowerInfo(ti)
  , _time(ti._time)
  , _amplitude(ti._amplitude)
{
}

TowerInfov1::~TowerInfov1()
{
}

void TowerInfov1::Reset()
{
  _time = 0;
  _amplitude = 0;
}
