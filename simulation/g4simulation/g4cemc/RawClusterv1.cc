#include "RawClusterv1.h"
#include <limits>

using namespace std;

RawClusterv1::RawClusterv1()
  : RawCluster()
  , clusterid(0)
  , _z(numeric_limits<float>::signaling_NaN())
  , _r(numeric_limits<float>::signaling_NaN())
  , _phi(numeric_limits<float>::signaling_NaN())
  , _energy(numeric_limits<float>::signaling_NaN())
{
}

void RawClusterv1::Reset()
{
  clusterid = 0;
  _z = (numeric_limits<float>::signaling_NaN());
  _r = (numeric_limits<float>::signaling_NaN());
  _phi = (numeric_limits<float>::signaling_NaN());
  _energy = (numeric_limits<float>::signaling_NaN());
  towermap.clear();
}

void RawClusterv1::addTower(const RawClusterDefs::keytype twrid, const float etower)
{
  if (towermap.find(twrid) != towermap.end())
  {
    cout << "tower 0x" << hex << twrid << ", dec: " << dec
         << twrid << " already exists, that is bad" << endl;
    exit(1);
  }
  towermap[twrid] = etower;
}

void RawClusterv1::identify(std::ostream& os = std::cout) const
{
  os << "RawClusterv1" << std::endl;
}

//! convert cluster location to psuedo-rapidity given a user chosen z-location
float RawClusterv1::get_eta(const float z) const
{
  if (get_r() <= 0) return numeric_limits<float>::signaling_NaN();
  return asinh((get_z() - z) / get_r());
}


//! convert cluster E_T given a user chosen z-location
float RawClusterv1::get_et(const float z) const
{
  if (get_r() <= 0) return numeric_limits<float>::signaling_NaN();
  return asinh((get_z() - z) / get_r());
}
