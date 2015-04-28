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


RawTowerv2* RawTowerv2::clone()
{
  int eta_i = this->get_bineta();
  int phi_i = this->get_binphi();

  RawTowerv2* tower_clone = new RawTowerv2( eta_i , phi_i );

  tower_clone->set_edep( this->get_edep() );
  tower_clone->set_thetaMin( this->get_thetaMin() );
  tower_clone->set_thetaSize( this->get_thetaSize() );
  tower_clone->set_phiMin( this->get_phiMin() );
  tower_clone->set_phiSize( this->get_phiSize() );
  tower_clone->set_zMin( this->get_zMin() );
  tower_clone->set_zSize( this->get_zSize() );

  return tower_clone;
}
