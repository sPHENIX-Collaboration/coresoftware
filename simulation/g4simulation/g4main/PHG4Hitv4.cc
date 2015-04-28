#include "PHG4Hitv4.h"
#include "PHG4HitDefs.h"

using namespace std;

G4Allocator<PHG4Hitv4> PHG4Hitv4Allocator;

ClassImp(PHG4Hitv4)

PHG4Hitv4::PHG4Hitv4():
 eion(NAN)
{}

PHG4Hitv4::PHG4Hitv4(PHG4Hit const &g4hit)
{
  Copy(g4hit);
}

void
PHG4Hitv4::print() const {
  std::cout<<"New Hitv4   "<<hitid<<" layer "<<layer<<"  on track "<<trackid<<" EDep "<<edep
	   <<" Eion: " << eion <<std::endl;
  std::cout<<"Location: X "<<x[0]<<"/"<<x[1]<<"  Y "<<y[0]<<"/"<<y[1]<<"  Z "<<z[0]<<"/"<<z[1]<<std::endl;
  std::cout<<"Time        "<<t[0]<<"/"<<t[1]<<std::endl;          
}
