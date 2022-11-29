#include "TowerInfo.h"

#include <phool/PHObject.h>  // for PHObject

#include <algorithm>
#include <cassert>
#include <cmath>

TowerInfo::TowerInfo()
{

}

TowerInfo::TowerInfo(const TowerInfo &ti)
  : TObject(ti)
  , _time(ti._time)
  , _amplitude(ti._amplitude)
{

}

TowerInfo::~TowerInfo()
{


}
