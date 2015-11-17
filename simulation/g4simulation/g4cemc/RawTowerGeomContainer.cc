#include "RawTowerGeomContainer.h"
#include "RawTowerGeom.h"

#include <cstdlib>
#include <iostream>

ClassImp(RawTowerGeomContainer)

using namespace std;



void
RawTowerGeomContainer::identify(std::ostream& os) const
{
  os << "Base class RawTowerGeomContainer."<< std::endl;
}
