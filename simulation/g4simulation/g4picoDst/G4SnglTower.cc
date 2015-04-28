#include "G4SnglTower.h"

using namespace std;

ClassImp(G4SnglTower)

G4SnglTower::G4SnglTower()
{
	Reset();
}

G4SnglTower::~G4SnglTower()
{
}

void G4SnglTower::Reset()
{
	detid = -9999;
	ieta = -9999.;
	iphi = -9999.;
	eta = -9999.;
	phi = -9999.;
	edep = -9999.;
}

void G4SnglTower::print() const
{
	std::cout<<"-----------------------------G4SnglTower::print()-----------------------------"<<endl;
	std::cout<<"detector id: "<<detid<<"   ieta: "<<ieta<<"  eta: "<<eta<<"  iphi: "<<iphi<<"  phi: "
	<<phi<<"  edep: "<<edep<<std::endl;
}
