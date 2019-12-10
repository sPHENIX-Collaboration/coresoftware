#include "PHG4Cellv2.h"

#include <g4main/PHG4HitDefs.h>  // for keytype

#include <iostream>

using namespace std;

PHG4Cellv2::PHG4Cellv2()
  : PHG4Cellv2(~0x0)
{
}

PHG4Cellv2::PHG4Cellv2(const PHG4CellDefs::keytype g4cellid)
  : cellid(g4cellid)
  , _edep(0)
{
}

PHG4Cellv2::~PHG4Cellv2()
{
  //  hitedeps.clear();
  return;
}

void PHG4Cellv2::add_edep(const PHG4HitDefs::keytype g4hitid, const float edep)
{
  hitedeps[g4hitid] = edep;
  return;
}

bool PHG4Cellv2::has_binning(const PHG4CellDefs::CellBinning binning) const
{
  return PHG4CellDefs::has_binning(cellid, binning);
}

short int
PHG4Cellv2::get_detid() const
{
  return PHG4CellDefs::get_detid(cellid);
}

void PHG4Cellv2::print() const
{
  identify(cout);
}

void PHG4Cellv2::Reset()
{
  hitedeps.clear();
  return;
}

void PHG4Cellv2::identify(std::ostream& os) const
{
  os << "New PHG4Cellv2  0x" << hex << cellid << dec << endl;

  os << "Associated to " << hitedeps.size() << " hits" << endl;
  for (const auto pair : hitedeps)
  {
    os << "\t PHG4hit " << pair.first << " -> " << pair.second << " GeV" << endl;
  }

  //  os <<"Contains to "<<trainOfDigits.size()<<" TPC digitization chain"<<endl;
}
