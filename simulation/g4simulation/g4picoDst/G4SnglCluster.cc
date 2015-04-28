#include "G4SnglCluster.h"

using namespace std;

ClassImp(G4SnglCluster)

G4SnglCluster::G4SnglCluster()
{
	Reset();
}

G4SnglCluster::~G4SnglCluster()
{
}

void G4SnglCluster::Reset()
{
	detid = -9999;
	ntowers = -9999;
	eta = -9999.;
	phi = -9999.;
	edep = -9999.;
}

void G4SnglCluster::print() const
{
	std::cout<<"-----------------------------G4SnglCluster::print()-----------------------------"<<endl;
	std::cout<<"detector id: "<<detid<<"    eta: "<<eta<<"    phi: "<<phi<<"    ntowers: "
	<< ntowers <<"    edep: "<<edep<<std::endl;
}
