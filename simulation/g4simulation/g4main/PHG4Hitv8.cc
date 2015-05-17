#include "PHG4Hitv8.h"

using namespace std;

G4Allocator<PHG4Hitv8> PHG4Hitv8Allocator;

ClassImp(PHG4Hitv8)

PHG4Hitv8::PHG4Hitv8()
{
  set_index_j(-1);
  set_index_k(-1);
  set_index_l(-1);
}

PHG4Hitv8::PHG4Hitv8(PHG4Hit const &g4hit)
{
  Copy(g4hit);
}

void
PHG4Hitv8::print() const {
  cout<<"New Hitv8    "<<hitid<<" layer "<<layer<<"  on track "<<trackid<<" EDep "<<edep<<endl;
  cout<<"Location: X  "<<x[0]<<"/"<<x[1]<<"  Y "<<y[0]<<"/"<<y[1]<<"  Z "<<z[0]<<"/"<<z[1]<<endl;
  cout<<"Time         "<<t[0]<<"/"<<t[1]<<endl;
  cout<<"Index j = "<< _idx_j <<" k = "<< _idx_k <<" l = "<< _idx_l <<endl;
}
