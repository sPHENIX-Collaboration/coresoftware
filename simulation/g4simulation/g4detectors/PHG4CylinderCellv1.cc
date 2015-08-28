#include "PHG4CylinderCellv1.h"
#include "PHG4CylinderCellDefs.h"

using namespace std;

ClassImp(PHG4CylinderCellv1)

PHG4CylinderCellv1::PHG4CylinderCellv1():
  layer(0xFFFFFFFF),
  cellid(0xFFFFFFFF),
  binz(-1),
  binphi(-1),
  edeps(), light_yield(0)
{}

void
PHG4CylinderCellv1::add_edep(const PHG4HitDefs::keytype g4hitid, const float edep)
{
  if (edeps.find(g4hitid) == edeps.end())
    {
      edeps[g4hitid] = edep;
    }
  else
    {
      edeps[g4hitid]+= edep;
    }
}

void
PHG4CylinderCellv1::add_edep(const PHG4HitDefs::keytype g4hitid, const float edep, const float ly)
{
  add_edep(g4hitid, edep);
  light_yield += ly;
}

double PHG4CylinderCellv1::get_edep() const {
  
  double esum = 0.0;
  EdepConstIterator iter;
  for (iter = edeps.begin(); iter != edeps.end(); ++iter) {
    esum += iter->second;
  }
  return esum;
}

void
PHG4CylinderCellv1::identify(std::ostream& os) const
{
  os << "PHG4CylinderCellv1: #" << cellid << " ";
  os << "(layer,binz,binphi,e) = (";
  os << layer << ","; 
  os << binz << ",";
  os << binphi << ",";
  os << get_edep();
  os << ")";
  os << endl;
}
