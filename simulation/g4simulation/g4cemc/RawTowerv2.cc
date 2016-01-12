#include "RawTowerv2.h"
#include <iostream>
#include <algorithm>

#include <map>

using namespace std;

ClassImp(RawTowerv2)

RawTowerv2::RawTowerv2() : bineta(-1), binphi(-1)
{
  thetaMin = -1;
  thetaSize = -1;
  phiMin = -1;
  phiSize = -1;
  zMin = -1;
  zSize = -1;
  edep = -1;
}

RawTowerv2::RawTowerv2(const int ieta, const int iphi) :
  bineta(ieta),
  binphi(iphi)
{
  thetaMin = -1;
  thetaSize = -1;
  phiMin = -1;
  phiSize = -1;
  zMin = -1;
  zSize = -1;
  edep = -1;
}

RawTowerv2::~RawTowerv2()
{}

void RawTowerv2::identify(std::ostream& os) const
{
  os << "RawTowerv2: etabin=(" << bineta << "), phibin=(" << binphi << "),energy=" << get_energy() << std::endl;

  os << "RawTowerv2: theta from " << thetaMin << " to " << thetaMin+thetaSize << "; phi from " << phiMin << " to " << phiMin+phiSize << "; energy " << get_energy() << std::endl;

}


RawTowerv2* RawTowerv2::clone() const
{
  // no dynamic memory allocation, just use the compiler provided
  // copy constructor to support this task
  // will prevent the "I forgot to add the new thing here" bugs
  RawTowerv2* tower_clone = new RawTowerv2(*this);
  return tower_clone;
}
