// $Id: $                                                                                             

/*!
 * \file PHG4CylinderCellSpacalv1.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4CylinderCell_Spacalv1.h"


ClassImp(PHG4CylinderCell_Spacalv1);

PHG4CylinderCell_Spacalv1::PHG4CylinderCell_Spacalv1():
    fiber_ID(-1)
{
  // TODO Auto-generated constructor stub

}

PHG4CylinderCell_Spacalv1::~PHG4CylinderCell_Spacalv1()
{
  // TODO Auto-generated destructor stub
}


pair<double, double>
PHG4CylinderCell_Spacalv1::get_zbounds(const int ibin) const
{
  if (ibin < 0 || ibin > nzbins)
    {
      cout << "Asking for invalid bin in z: " << ibin << endl;
      exit(1);
    }
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  double zlow = zmin + ibin * zstep;
  double zhigh = zlow + zstep;
  return make_pair(zlow, zhigh);
}

pair<double, double>
PHG4CylinderCell_Spacalv1::get_etabounds(const int ibin) const
{
  if (ibin < 0 || ibin > nzbins)
    {
      cout << "Asking for invalid bin in z: " << ibin << endl;
      exit(1);
    }
  check_binning_method_eta("PHG4CylinderCell_Spacalv1::get_etabounds");
//  check_binning_method(PHG4CylinderCellDefs::etaphibinning);
  double zlow = zmin + ibin * zstep;
  double zhigh = zlow + zstep;
  return make_pair(zlow, zhigh);
}


pair<double, double>
PHG4CylinderCell_Spacalv1::get_phibounds(const int ibin) const
{
  if (ibin < 0 || ibin > nphibins)
    {
      cout << "Asking for invalid bin in phi: " << ibin << endl;
      exit(1);
    }

  double philow = phimin + ibin * phistep;
  double phihigh = philow + phistep;
  return make_pair(philow, phihigh);
}

int
PHG4CylinderCell_Spacalv1::get_zbin(const double z) const
{
  if (z < zmin || z > (zmin+nzbins*zstep))
  {
    //    cout << "Asking for bin for z outside of z range: " << z << endl;
    return -1;
  }

  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  return floor( (z-zmin)/zstep );
}

int
PHG4CylinderCell_Spacalv1::get_etabin(const double eta) const
{
  if (eta < zmin || eta > (zmin+nzbins*zstep))
  {
    //    cout << "Asking for bin for eta outside of eta range: " << eta << endl;
    return -1;
  }
  check_binning_method_eta();
  return floor( (eta-zmin)/zstep );
}

int
PHG4CylinderCell_Spacalv1::get_phibin(const double phi) const
{
  double norm_phi = phi;
  if(phi < phimin || phi > (phimin+nphibins*phistep))
  {
    int nwraparound = -floor((phi-phimin) * 0.5 / M_PI);
    norm_phi += 2*M_PI*nwraparound;
  }
  check_binning_method_phi();
  return floor( (norm_phi-phimin)/phistep );
}

double
PHG4CylinderCell_Spacalv1::get_zcenter(const int ibin) const
{
  if (ibin < 0 || ibin > nzbins)
    {
      cout << "Asking for invalid bin in z: " << ibin << endl;
      exit(1);
    }
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  return zmin + (ibin + 0.5)*zstep;
}

double
PHG4CylinderCell_Spacalv1::get_etacenter(const int ibin) const
{
  if (ibin < 0 || ibin > nzbins)
    {
      cout << "Asking for invalid bin in eta: " << ibin << endl;
      cout << "minbin: 0, maxbin " << nzbins << endl;
      exit(1);
    }
  check_binning_method_eta();
  return zmin + (ibin + 0.5)*zstep;
}

double
PHG4CylinderCell_Spacalv1::get_phicenter(const int ibin) const
{
  if (ibin < 0 || ibin > nphibins)
    {
      cout << "Asking for invalid bin in phi: " << ibin << endl;
      exit(1);
    }

  check_binning_method_phi();
  return (phimin + (ibin + 0.5)*phistep);
}
