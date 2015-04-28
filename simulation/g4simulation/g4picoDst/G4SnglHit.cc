#include "G4SnglHit.h"

using namespace std;

ClassImp(G4SnglHit)

G4SnglHit::G4SnglHit()
{
	Reset();
}

G4SnglHit::~G4SnglHit()
{
}

void G4SnglHit::Reset()
{
	detid = -9999;
	layer = -9999;
	scintid = -9999;
	trackid = -9999;
	x0 = -9999.;
	y0 = -9999.;
	z0 = -9999.;
	x1 = -9999.;
	y1 = -9999.;
	z1 = -9999.;
	edep = -9999.;
}

void G4SnglHit::print() const
{
	std::cout<<"-----------------------------G4SnglHit::print()-----------------------------"<<endl;
	std::cout<<"detector id: "<<detid<<"   scint id: "<<scintid<<"   layer: "<<layer<<"   track id: "<<trackid<<"   EDep: "<<edep<<std::endl;
	std::cout<<"position: X "<<x0<<"/"<<x1<<"  Y "<<y0<<"/"<<y1<<"  Z "<<z0<<"/"<<z1<<std::endl;
}
