// $Id: $                                                                                             

/*!
 * \file PHG4CylinderCellGeomSpacalv1.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4CylinderCellGeom_Spacalv1.h"

ClassImp(PHG4CylinderCellGeom_Spacalv1);

using namespace std;

PHG4CylinderCellGeom_Spacalv1::PHG4CylinderCellGeom_Spacalv1()
{
  // TODO Auto-generated constructor stub

}

PHG4CylinderCellGeom_Spacalv1::~PHG4CylinderCellGeom_Spacalv1()
{
  // TODO Auto-generated destructor stub
}

void
PHG4CylinderCellGeom_Spacalv1::map_consistency_check()
{
  if (nzbins != z_bound_map.size())
    {
      cout << "PHG4CylinderCellGeom_Spacalv1::map_consistency_check - "
          << "z_bound_map.size() of " << z_bound_map.size()
          << " in inconsistent with nzbins of " << nzbins << endl;
      exit(1);
    }
  if (nzbins != eta_bound_map.size())
    {
      cout << "PHG4CylinderCellGeom_Spacalv1::map_consistency_check - "
          << "eta_bound_map.size() of " << eta_bound_map.size()
          << " in inconsistent with nzbins of " << nzbins << endl;
      exit(1);
    }
}

pair<double, double>
PHG4CylinderCellGeom_Spacalv1::get_zbounds(const int ibin) const
{
  if (ibin < 0 || ibin > nzbins)
    {
      cout
          << "PHG4CylinderCellGeom_Spacalv1::get_zbounds - Asking for invalid bin in z: "
          << ibin << endl;
      exit(1);
    }
  map_consistency_check();
  check_binning_method(PHG4CylinderCellDefs::spacalbinning);
  return z_bound_map[ibin];
}

pair<double, double>
PHG4CylinderCellGeom_Spacalv1::get_etabounds(const int ibin) const
{
  if (ibin < 0 || ibin > nzbins)
    {
      cout
          << "PHG4CylinderCellGeom_Spacalv1::get_etabounds - Asking for invalid bin in z: "
          << ibin << endl;
      exit(1);
    }
  map_consistency_check();
  check_binning_method(PHG4CylinderCellDefs::spacalbinning);
  return eta_bound_map[ibin];
}

int
PHG4CylinderCellGeom_Spacalv1::get_zbin(const double z) const
{
  cout << "PHG4CylinderCellGeom_Spacalv1::get_zbin is invalid" << endl;
  exit(1);
  return -1;
}

int
PHG4CylinderCellGeom_Spacalv1::get_etabin(const double eta) const
{
  cout << "PHG4CylinderCellGeom_Spacalv1::get_etabin is invalid" << endl;
  exit(1);
  return -1;
}

double
PHG4CylinderCellGeom_Spacalv1::get_zcenter(const int ibin) const
{
  if (ibin < 0 || ibin > nzbins)
    {
      cout
          << "PHG4CylinderCellGeom_Spacalv1::get_zcenter - Asking for invalid bin in z: "
          << ibin << endl;
      exit(1);
    }
  map_consistency_check();
  check_binning_method(PHG4CylinderCellDefs::spacalbinning);
  return 0.5 * (z_bound_map[ibin].first + z_bound_map[ibin].second);
}

double
PHG4CylinderCellGeom_Spacalv1::get_etacenter(const int ibin) const
{
  if (ibin < 0 || ibin > nzbins)
    {
      cout << "Asking for invalid bin in eta: " << ibin << endl;
      cout << "minbin: 0, maxbin " << nzbins << endl;
      exit(1);
    }
  return 0.5 * (eta_bound_map[ibin].first + eta_bound_map[ibin].second);
}
