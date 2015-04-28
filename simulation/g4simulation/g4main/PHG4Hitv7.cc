#include "PHG4Hitv7.h"
#include "PHG4HitDefs.h"

using namespace std;

G4Allocator<PHG4Hitv7> PHG4Hitv7Allocator;

ClassImp(PHG4Hitv7)

PHG4Hitv7::PHG4Hitv7()
{
  set_strip_z_index(-9999);
  set_strip_y_index(-9999);
  set_ladder_z_index(-9999);
  set_ladder_phi_index(-9999);
}

PHG4Hitv7::PHG4Hitv7(PHG4Hit const &g4hit)
{
  Copy(g4hit);
}

void
PHG4Hitv7::print() const {
  cout<<"New Hitv7    "<<hitid<<" layer "<<layer<<"  on track "<<trackid<<" EDep "<<edep<<endl;
  cout<<"Location: X  "<<x[0]<<"/"<<x[1]<<"  Y "<<y[0]<<"/"<<y[1]<<"  Z "<<z[0]<<"/"<<z[1]<<endl;
  cout<<"Time         "<<t[0]<<"/"<<t[1]<<endl; 
  cout<<"Momentum  Px "<<px[0]<<"/"<<px[1]<<"  Py "<<py[0]<<"/"<<py[1]<<"  Pz "<<pz[0]<<"/"<<pz[1]<<endl; 
  cout<<"strip_z_index "<<strip_z_index<<" strip_y_index"<<strip_y_index<<"  ladder_z_index "<<ladder_z_index<<"  ladder_phi_index "<<ladder_phi_index<<endl; 
}

