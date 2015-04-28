#include "PHG4Hitv3.h"

using namespace std;

G4Allocator<PHG4Hitv3> PHG4Hitv3Allocator;

ClassImp(PHG4Hitv3)

  PHG4Hitv3::PHG4Hitv3():
    scint_id(-9999)
{}

PHG4Hitv3::PHG4Hitv3(PHG4Hit const &g4hit)
{
  Copy(g4hit);
}

void
PHG4Hitv3::print() const {
  cout<<"New Hitv3    "<<hitid<<" layer "<<layer<<"  on track "<<trackid<<" EDep "<<edep<<endl;
  cout<<"Location: X  "<<x[0]<<"/"<<x[1]<<"  Y "<<y[0]<<"/"<<y[1]<<"  Z "<<z[0]<<"/"<<z[1]<<endl;
  cout<<"Time         "<<t[0]<<"/"<<t[1]<<endl; 
  cout<<"Scintillator id "<< scint_id << endl;
}

