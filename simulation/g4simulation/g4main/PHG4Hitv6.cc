// $Id: PHG4Hitv6.C,v 1.1 2015/01/06 02:33:45 jinhuang Exp $                                                                                             

/*!
 * \file PHG4Hitv6.C
 * \brief PHG4Hitv6 with scintillation light yield and pathlength
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.1 $
 * \date $Date: 2015/01/06 02:33:45 $
 */

#include "PHG4Hitv6.h"

using namespace std;

G4Allocator<PHG4Hitv6> PHG4Hitv6Allocator;

ClassImp(PHG4Hitv6)

PHG4Hitv6::PHG4Hitv6() :
    light_yield(NAN), path_length(NAN)
{
}

PHG4Hitv6::~PHG4Hitv6()
{
  // TODO Auto-generated destructor stub
}

PHG4Hitv6::PHG4Hitv6(PHG4Hit const &g4hit)
{
  Copy(g4hit);
}

void
PHG4Hitv6::print() const
{
  cout << "New Hitv5    " << hitid << " layer " << layer << "  on track "
      << trackid << " EDep " << edep << endl;
  cout << "Location: X  " << x[0] << "/" << x[1] << "  Y " << y[0] << "/"
      << y[1] << "  Z " << z[0] << "/" << z[1] << endl;
  cout << "Time         " << t[0] << "/" << t[1] << endl;
  cout << "Path Length  " << path_length << endl;
  cout << "Light Yield  " << light_yield << endl;
  cout << "Scintillator id " << scint_id << endl;
}
