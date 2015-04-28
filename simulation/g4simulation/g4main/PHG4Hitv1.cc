#include "PHG4Hitv1.h"
#include "PHG4HitDefs.h"

using namespace std;

G4Allocator<PHG4Hitv1> PHG4Hitv1Allocator;

ClassImp(PHG4Hitv1)

PHG4Hitv1::PHG4Hitv1():
 layer(0xffffffff),
 hitid(0xffffffff),
 trackid(-9999),
 edep(NAN)
{
  for (int i = 0; i<2;i++)
    {
      set_x(i,NAN);
      set_y(i,NAN);
      set_z(i,NAN);
      set_t(i,NAN);
    }
}

PHG4Hitv1::PHG4Hitv1(PHG4Hit const &g4hit)
{
  Copy(g4hit);
}

void
PHG4Hitv1::print() const {
  std::cout<<"New Hitv1   "<<hitid<<" layer "<<layer<<"  on track "<<trackid<<" EDep "<<edep<<std::endl;
  std::cout<<"Location: X "<<x[0]<<"/"<<x[1]<<"  Y "<<y[0]<<"/"<<y[1]<<"  Z "<<z[0]<<"/"<<z[1]<<std::endl;
  std::cout<<"Time        "<<t[0]<<"/"<<t[1]<<std::endl;          
}

