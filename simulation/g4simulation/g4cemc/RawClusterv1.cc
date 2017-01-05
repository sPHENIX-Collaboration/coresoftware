#include "RawClusterv1.h"

using namespace std;

RawClusterv1::RawClusterv1():
  RawCluster(),
  clusterid(0),
  _eta(0.0),
  _phi(0.0),
  _energy(0.0)
{}

void
RawClusterv1::Reset()
{
  clusterid = 0;
  _eta = 0.0;
  _phi = 0.0;
  _energy = 0.0;
  towermap.clear();
}

void
RawClusterv1::addTower(const RawClusterDefs::keytype twrid, const float etower)
{
  if (towermap.find(twrid) != towermap.end())
    {
      cout << "tower 0x" << hex << twrid << ", dec: " << dec
	   << twrid << " already exists, that is bad" << endl;
      exit(1);
    }
  towermap[twrid] = etower;
}
