#include "RawTowerGeomv2.h"

#include <cmath>
#include <cstdlib>
#include <cassert>

ClassImp(RawTowerGeomv2)

using namespace std;

RawTowerGeomv2::RawTowerGeomv2() :
    radius(NAN), thickness(NAN), nphibins(-1), phimin(-M_PI), phistep(NAN)
{
  return;
}

RawTowerGeomv2::RawTowerGeomv2(RawTowerGeom *geo)
{
  geo->identify();
  set_radius(geo->get_radius());
  set_thickness(geo->get_thickness());
  set_phibins(geo->get_phibins());
  set_phistep(geo->get_phistep());
  set_phimin(geo->get_phimin());
  set_etabins(geo->get_etabins());
  for (int i=0; i<geo->get_etabins(); i++)
    {
      set_etabounds(i,geo->get_etabounds(i));
    }
}

void
RawTowerGeomv2::set_etabins(const int i)
{
  assert(i>0);
  bound_t invalid_bound(NAN, NAN);
  eta_bound_map.resize(i, invalid_bound);
}

void
RawTowerGeomv2::identify(std::ostream& os) const
{
  os << "RawTowerGeomv2: radius: " << radius << ", thickness: " << thickness
      << ", etabins: " << get_etabins() << ", phimin: " << phimin
      << ", phibins: " << nphibins << ", phistep: " << phistep
      << ", eta bins are: ";

  int i = 0;
  for (bound_map_t::const_iterator iter = eta_bound_map.begin();
      iter != eta_bound_map.end(); ++iter)
    {
      cout << "bin[" << i << "](" << iter->first << ", " << iter->second
          << ")  ";
      i++;
    }
  cout << endl;

  return;
}

pair<double, double>
RawTowerGeomv2::get_etabounds(const int ibin) const
{
  if (ibin < 0 || ibin > get_etabins())
    {
      cout << "Asking for invalid bin in eta: " << ibin << endl;
      exit(1);
    }
  return eta_bound_map[ibin];
}

pair<double, double>
RawTowerGeomv2::get_phibounds(const int ibin) const
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
RawTowerGeomv2::get_etabin(const double eta) const
{

  int ibin = -1;
  int i = 0;
  for (bound_map_t::const_iterator iter = eta_bound_map.begin();
      iter != eta_bound_map.end(); ++iter)
    {
      if (eta >= iter->first && eta < iter->second)
        {
          ibin = i;
          break;
        }

      i++;
    }
  if (ibin < 0 || ibin >= nphibins)
    {
      cout
          << "RawTowerGeomv2::get_etabin - ERROR - Asking for invalid bin in eta "
          << eta << endl;
      exit(1);
    }

  return ibin;
}

int
RawTowerGeomv2::get_phibin(const double phi) const
{
  double norm_phi = phi;
  if (phi < phimin || phi > (phimin + nphibins * phistep))
    {
      int nwraparound = -floor((phi - phimin) * 0.5 / M_PI);
      norm_phi += 2 * M_PI * nwraparound;
    }
  return floor((norm_phi - phimin) / phistep);
}

double
RawTowerGeomv2::get_etacenter(const int ibin) const
{
  if (ibin < 0 || ibin >= get_etabins())
    {
      cout << "RawTowerGeomv2::get_etacenter - Asking for invalid bin in eta: "
	   << ibin << endl;
      cout << "minbin: 0, maxbin " << get_etabins() << endl;
      exit(1);
    }
  return (eta_bound_map[ibin].first + eta_bound_map[ibin].second) / 2.;
}

double
RawTowerGeomv2::get_phicenter(const int ibin) const
{
  if (ibin < 0 || ibin > nphibins)
    {
      cout << "Asking for invalid bin in phi: " << ibin << endl;
      exit(1);
    }

  return (phimin + (ibin + 0.5) * phistep);
}

void
RawTowerGeomv2::
set_etabounds(const int ibin, const std::pair<double, double> & bounds)
{

  if (ibin < 0 || ibin >= get_etabins())
    {
      cout << "RawTowerGeomv2::set_bounds - Asking for invalid bin in eta: "
	   << ibin << endl;
      cout << "minbin: 0, maxbin " << get_etabins() << endl;
      exit(1);
    }

  std::pair<double, double> b_reg(bounds);
  if (b_reg.first>b_reg.second)
    {
      b_reg.second = bounds.first;
      b_reg.first = bounds.second;
    }

  eta_bound_map[ibin] = b_reg;
}
