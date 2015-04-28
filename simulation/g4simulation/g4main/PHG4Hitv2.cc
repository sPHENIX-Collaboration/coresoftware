#include "PHG4Hitv2.h"
#include "PHG4HitDefs.h"

using namespace std;

G4Allocator<PHG4Hitv2> PHG4Hitv2Allocator;

ClassImp(PHG4Hitv2)

PHG4Hitv2::PHG4Hitv2()
{
  for (int i = 0; i<2;i++)
    {
      set_px(i,NAN);
      set_py(i,NAN);
      set_pz(i,NAN);
    }
}

PHG4Hitv2::PHG4Hitv2(PHG4Hit const &g4hit)
{
  Copy(g4hit);
}

void
PHG4Hitv2::print() const {
  cout<<"New Hitv2    "<<hitid<<" layer "<<layer<<"  on track "<<trackid<<" EDep "<<edep<<endl;
  cout<<"Location: X  "<<x[0]<<"/"<<x[1]<<"  Y "<<y[0]<<"/"<<y[1]<<"  Z "<<z[0]<<"/"<<z[1]<<endl;
  cout<<"Time         "<<t[0]<<"/"<<t[1]<<endl; 
  cout<<"Momentum  Px "<<px[0]<<"/"<<px[1]<<"  Py "<<py[0]<<"/"<<py[1]<<"  Pz "<<pz[0]<<"/"<<pz[1]<<endl; 
}

