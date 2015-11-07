#include "RawTowerGeomv1.h"

#include <cmath>
#include <cstdlib>

ClassImp(RawTowerGeomv1)

using namespace std;

RawTowerGeomv1::RawTowerGeomv1():
  radius(NAN),
  thickness(NAN),
  netabins(-1),
  etamin(NAN),
  etastep(NAN),
  nphibins(-1),
  phimin(-M_PI),
  phistep(NAN)
{
  return;
}

void
RawTowerGeomv1::identify(std::ostream& os) const
{
  os << "RawTowerGeomv1: radius: " << radius
     << ", thickness: " << thickness
     << ", etabins: " << netabins
     << ", etamin: " << etamin
     << ", etastepsize: " << etastep
     << ", phimin: " << phimin
     << ", phibins: " << nphibins
     << ", phistep: " << phistep
     << endl;
  return;
}

pair<double, double>
RawTowerGeomv1::get_etabounds(const int ibin) const
{
  if (ibin < 0 || ibin > netabins)
    {
      cout << "Asking for invalid bin in eta: " << ibin << endl;
      exit(1);
    }
  double etalow = etamin + ibin * etastep;
  double etahigh = etalow + etastep;
  return make_pair(etalow, etahigh);
}


pair<double, double>
RawTowerGeomv1::get_phibounds(const int ibin) const
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
RawTowerGeomv1::get_etabin(const double eta) const
{
  if (eta < etamin || eta > (etamin+netabins*etastep))
  {
    cout << "Asking for bin for eta outside of eta range: " << eta << endl;
    return -1;
  }
  return floor( (eta-etamin)/etastep );
}

int
RawTowerGeomv1::get_phibin(const double phi) const
{
//  double norm_phi = phi;
//  if (phi < phimin || phi > (phimin + nphibins * phistep))
//    {
//      int nwraparound = -floor((phi - phimin) * 0.5 / M_PI);
//      norm_phi += 2 * M_PI * nwraparound;
//    }

  const int bin = floor((phi - phimin) / phistep);
  const int bin_wrap = floor((double)bin/(double)nphibins)*nphibins;

  return bin - bin_wrap;
}


double
RawTowerGeomv1::get_etacenter(const int ibin) const
{
  if (ibin < 0 || ibin > netabins)
    {
      cout << "Asking for invalid bin in eta: " << ibin << endl;
      cout << "minbin: 0, maxbin " << netabins << endl; 
      exit(1);
    }
  return etamin + (ibin + 0.5)*etastep;
}

double
RawTowerGeomv1::get_phicenter(const int ibin) const
{
  if (ibin < 0 || ibin > nphibins)
    {
      cout << "Asking for invalid bin in phi: " << ibin << endl;
      exit(1);
    }

  return (phimin + (ibin + 0.5)*phistep);
}

