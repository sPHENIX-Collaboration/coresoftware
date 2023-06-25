#include "RawTowerGeomContainer_Cylinderv1.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>

using namespace std;

RawTowerGeomContainer_Cylinderv1::RawTowerGeomContainer_Cylinderv1(RawTowerDefs::CalorimeterId caloid)
  : RawTowerGeomContainerv1(caloid)
  , radius(NAN)
  , thickness(NAN)
{
  return;
}

void RawTowerGeomContainer_Cylinderv1::Reset()
{
  eta_bound_map.clear();

  phi_bound_map.clear();

  RawTowerGeomContainerv1::Reset();
}

void RawTowerGeomContainer_Cylinderv1::set_etabins(const int i)
{
  assert(i > 0);
  bound_t invalid_bound(NAN, NAN);
  eta_bound_map.resize(i, invalid_bound);
}

void RawTowerGeomContainer_Cylinderv1::set_phibins(const int i)
{
  assert(i > 0);
  bound_t invalid_bound(NAN, NAN);
  phi_bound_map.resize(i, invalid_bound);
}

void RawTowerGeomContainer_Cylinderv1::identify(std::ostream& os) const
{
  os << "RawTowerGeomContainer_Cylinderv1: radius: " << radius
     << ", thickness: " << thickness << ", etabins: " << get_etabins()
     << ", phibins: " << get_phibins();

  int i = 0;
  for (const auto& iter : eta_bound_map)
  {
    os << "eta_bin[" << i << "](" << iter.first << ", " << iter.second
       << ")  ";
    i++;
  }
  os << endl;
  i = 0;
  for (const auto& iter : phi_bound_map)
  {
    os << "phi_bin[" << i << "](" << iter.first << ", " << iter.second
       << ")  ";
    i++;
  }
  os << endl;
  return;
}

pair<double, double>
RawTowerGeomContainer_Cylinderv1::get_etabounds(const int ibin) const
{
  if (ibin < 0 || ibin > get_etabins())
  {
    identify();
    cout
        << "RawTowerGeomContainer_Cylinderv1::get_etabounds - Asking for invalid bin in eta: "
        << ibin << endl;
    exit(1);
  }
  return eta_bound_map[ibin];
}

pair<double, double>
RawTowerGeomContainer_Cylinderv1::get_phibounds(const int ibin) const
{
  if (ibin < 0 || ibin > get_phibins())
  {
    identify();
    cout
        << "RawTowerGeomContainer_Cylinderv1::get_phibounds - Asking for invalid bin in phi: "
        << ibin << endl;
    exit(1);
  }
  return phi_bound_map[ibin];
}

int RawTowerGeomContainer_Cylinderv1::get_etabin(const double eta) const
{
  int ibin = -1;
  int i = 0;

  // switch to search for the closest bin
  // since in a realistic calorimeter, there could be gaps
  double min_deta = 10;

  for (const auto& iter : eta_bound_map)
  {
    const double mean_eta = 0.5 * (iter.first + iter.second);

    if (eta >= iter.first && eta < iter.second)
    {
      // found the bin that the hit belong
      min_deta = 0;
      ibin = i;
      break;
    }
    else
    {
      const double deta = fabs(mean_eta - eta);
      if (deta < min_deta)
      {
        min_deta = deta;
        ibin = i;
      }  // keep searching
    }

    i++;
  }

  if (ibin < 0)
  {
    cout
        << "RawTowerGeomContainer_Cylinderv1::get_etabin - ERROR - Asking for invalid bin in eta "
        << eta << endl;
    exit(1);
  }

  return ibin;
}

int RawTowerGeomContainer_Cylinderv1::get_phibin(const double phi) const
{
  int ibin = -1;
  int i = 0;

  // switch to search for the closest bin
  // since in a realistic calorimeter, there could be gaps
  double min_dphi = 10;

  for (const auto& iter : phi_bound_map)
  {
    const double mean_phi = 0.5 * (iter.first + iter.second);

    const double phi_fold = phi - round((phi - mean_phi) / 2. / M_PI) * 2 * M_PI;

    if (phi_fold >= iter.first && phi_fold < iter.second)
    {
      // found the bin that the hit belong
      min_dphi = 0;
      ibin = i;
      break;
    }
    else
    {
      const double dphi = fabs(mean_phi - phi_fold);
      if (dphi < min_dphi)
      {
        min_dphi = dphi;
        ibin = i;
      }  // keep searching
    }

    i++;
  }

  if (ibin < 0)
  {
    cout
        << "RawTowerGeomContainer_Cylinderv1::get_phibin - ERROR - Asking for invalid bin in phi "
        << phi << endl;
    exit(1);
  }

  return ibin;
}

double
RawTowerGeomContainer_Cylinderv1::get_etacenter(const int ibin) const
{
  if (ibin < 0 || ibin >= get_etabins())
  {
    cout
        << "RawTowerGeomContainer_Cylinderv1::get_etacenter - Asking for invalid bin in eta: "
        << ibin << endl;
    cout << "minbin: 0, maxbin " << get_etabins() << endl;
    exit(1);
  }
  return (eta_bound_map[ibin].first + eta_bound_map[ibin].second) / 2.;
}

void RawTowerGeomContainer_Cylinderv1::set_etabounds(const int ibin,
                                                     const std::pair<double, double>& bounds)
{
  if (ibin < 0 || ibin >= get_etabins())
  {
    cout
        << "RawTowerGeomContainer_Cylinderv1::set_bounds - Asking for invalid bin in eta: "
        << ibin << endl;
    cout << "minbin: 0, maxbin " << get_etabins() << endl;
    exit(1);
  }

  std::pair<double, double> b_reg(bounds);
  if (b_reg.first > b_reg.second)
  {
    b_reg.second = bounds.first;
    b_reg.first = bounds.second;
  }

  eta_bound_map[ibin] = b_reg;
}

double
RawTowerGeomContainer_Cylinderv1::get_phicenter(const int ibin) const
{
  if (ibin < 0 || ibin >= get_phibins())
  {
    cout
        << "RawTowerGeomContainer_Cylinderv1::get_phicenter - Asking for invalid bin in phi: "
        << ibin << endl;
    cout << "minbin: 0, maxbin " << get_phibins() << endl;
    exit(1);
  }
  return (phi_bound_map[ibin].first + phi_bound_map[ibin].second) / 2.;
}

void RawTowerGeomContainer_Cylinderv1::set_phibounds(const int ibin,
                                                     const std::pair<double, double>& bounds)
{
  if (ibin < 0 || ibin >= get_phibins())
  {
    cout
        << "RawTowerGeomContainer_Cylinderv1::set_bounds - Asking for invalid bin in phi: "
        << ibin << endl;
    cout << "minbin: 0, maxbin " << get_phibins() << endl;
    exit(1);
  }

  std::pair<double, double> b_reg(bounds);
  if (b_reg.first > b_reg.second)
  {
    b_reg.second = bounds.first;
    b_reg.first = bounds.second;
  }

  phi_bound_map[ibin] = b_reg;
}
