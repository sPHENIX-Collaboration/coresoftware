#include "PHG4BlockCellGeom.h"
#include "PHG4CylinderCellDefs.h"

#include <cmath>
#include <cstdlib>

using namespace std;

PHG4BlockCellGeom::PHG4BlockCellGeom()
  : _layer(-9999)
  , _binning(0)
  , _nzbins(-1)
  , _zmin(NAN)
  , _zstep(NAN)
  , _nxbins(-1)
  , _xmin(NAN)
  , _xstep(NAN)
{
  return;
}

void PHG4BlockCellGeom::set_zbins(const int i)
{
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  _nzbins = i;
}

void PHG4BlockCellGeom::set_zmin(const double z)
{
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  _zmin = z;
}

int PHG4BlockCellGeom::get_zbins() const
{
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  return _nzbins;
}

double
PHG4BlockCellGeom::get_zmin() const
{
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  return _zmin;
}

double
PHG4BlockCellGeom::get_zstep() const
{
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  return _zstep;
}

void PHG4BlockCellGeom::set_zstep(const double z)
{
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  _zstep = z;
}

int PHG4BlockCellGeom::get_xbins() const
{
  check_binning_method_x("PHG4BlockCellGeom::get_xbins");
  return _nxbins;
}

double
PHG4BlockCellGeom::get_xstep() const
{
  check_binning_method_x("PHG4BlockCellGeom::get_xstep");
  return _xstep;
}

double
PHG4BlockCellGeom::get_xmin() const
{
  check_binning_method_x("PHG4BlockCellGeom::get_xmin");
  return _xmin;
}

void PHG4BlockCellGeom::set_xbins(const int i)
{
  check_binning_method_x("PHG4BlockCellGeom::set_xbins");
  _nxbins = i;
}

void PHG4BlockCellGeom::set_xstep(const double r)
{
  check_binning_method_x("PHG4BlockCellGeom::set_xstep");
  _xstep = r;
}

void PHG4BlockCellGeom::set_xmin(const double r)
{
  check_binning_method_x("PHG4BlockCellGeom::set_xmin");
  _xmin = r;
}

int PHG4BlockCellGeom::get_etabins() const
{
  check_binning_method_eta("PHG4BlockCellGeom::get_etabins");
  return _nzbins;
}

double
PHG4BlockCellGeom::get_etastep() const
{
  check_binning_method_eta("PHG4BlockCellGeom::get_etastep");
  return _zstep;
}
double
PHG4BlockCellGeom::get_etamin() const
{
  check_binning_method_eta("PHG4BlockCellGeom::get_etamin");
  return _zmin;
}

void PHG4BlockCellGeom::set_etamin(const double z)
{
  check_binning_method_eta("PHG4BlockCellGeom::set_etamin");
  _zmin = z;
}

void PHG4BlockCellGeom::set_etastep(const double z)
{
  check_binning_method_eta("PHG4BlockCellGeom::set_etastep");
  _zstep = z;
}

void PHG4BlockCellGeom::set_etabins(const int i)
{
  check_binning_method_eta("PHG4BlockCellGeom::set_etabins");
  _nzbins = i;
}

void PHG4BlockCellGeom::identify(std::ostream& os) const
{
  os << "layer: " << _layer;
  switch (_binning)
  {
  case PHG4CylinderCellDefs::sizebinning:
    os << ", zbins: " << _nzbins
       << ", zmin: " << _zmin
       << ", zstepsize: " << _zstep;
    break;

  case PHG4CylinderCellDefs::etaphibinning:
    os << ", etabins: " << _nzbins
       << ", etamin: " << _zmin
       << ", etastepsize: " << _zstep;
    break;

  case PHG4CylinderCellDefs::etaslatbinning:
    os << ", etabins: " << _nzbins
       << ", etamin: " << _zmin
       << ", etastepsize: " << _zstep;
    break;

  default:
    os << "no valid binning method: " << _binning << endl;
    return;
    break;
  }

  os << ", xmin: " << _xmin
     << ", xbins: " << _nxbins
     << ", xstep: " << _xstep
     << endl;
  return;
}

pair<double, double>
PHG4BlockCellGeom::get_zbounds(const int ibin) const
{
  if (ibin < 0 || ibin > _nzbins)
  {
    cout << "Asking for invalid bin in z: " << ibin << endl;
    exit(1);
  }
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  double zlow = _zmin + ibin * _zstep;
  double zhigh = zlow + _zstep;
  return make_pair(zlow, zhigh);
}

pair<double, double>
PHG4BlockCellGeom::get_etabounds(const int ibin) const
{
  if (ibin < 0 || ibin > _nzbins)
  {
    cout << "Asking for invalid bin in z: " << ibin << endl;
    exit(1);
  }
  check_binning_method_eta("PHG4BlockCellGeom::get_etabounds");
  //  check_binning_method(PHG4CylinderCellDefs::etaphibinning);
  double zlow = _zmin + ibin * _zstep;
  double zhigh = zlow + _zstep;
  return make_pair(zlow, zhigh);
}

pair<double, double>
PHG4BlockCellGeom::get_xbounds(const int ibin) const
{
  if (ibin < 0 || ibin > _nxbins)
  {
    cout << "Asking for invalid bin in x: " << ibin << endl;
    exit(1);
  }

  double xlow = _xmin + ibin * _xstep;
  double xhigh = xlow + _xstep;
  return make_pair(xlow, xhigh);
}

int PHG4BlockCellGeom::get_zbin(const double z) const
{
  if (z < _zmin || z > (_zmin + _nzbins * _zstep))
  {
    cout << "Asking for bin for z outside of z range: " << z << endl;
    return -1;
  }

  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  return floor((z - _zmin) / _zstep);
}

int PHG4BlockCellGeom::get_etabin(const double eta) const
{
  if (eta < _zmin || eta > (_zmin + _nzbins * _zstep))
  {
    cout << "Asking for bin for eta outside of eta range: " << eta << endl;
    return -1;
  }
  check_binning_method_eta();
  return floor((eta - _zmin) / _zstep);
}

int PHG4BlockCellGeom::get_xbin(const double x) const
{
  double norm_x = x;
  if (x < _xmin || x > (_xmin + _nxbins * _xstep))
  {
    cout << "Asking for bin for x outside of x range: " << x << endl;
    return -1;
  }
  check_binning_method_x();
  return floor((norm_x - _xmin) / _xstep);
}

double
PHG4BlockCellGeom::get_zcenter(const int ibin) const
{
  if (ibin < 0 || ibin > _nzbins)
  {
    cout << "Asking for invalid bin in z: " << ibin << endl;
    exit(1);
  }
  check_binning_method(PHG4CylinderCellDefs::sizebinning);
  return _zmin + (ibin + 0.5) * _zstep;
}

double
PHG4BlockCellGeom::get_etacenter(const int ibin) const
{
  if (ibin < 0 || ibin > _nzbins)
  {
    cout << "Asking for invalid bin in eta: " << ibin << endl;
    cout << "minbin: 0, maxbin " << _nzbins << endl;
    exit(1);
  }
  check_binning_method_eta();
  return _zmin + (ibin + 0.5) * _zstep;
}

double
PHG4BlockCellGeom::get_xcenter(const int ibin) const
{
  if (ibin < 0 || ibin > _nxbins)
  {
    cout << "Asking for invalid bin in x: " << ibin << endl;
    exit(1);
  }

  check_binning_method_x();
  return (_xmin + (ibin + 0.5) * _xstep);
}

string
PHG4BlockCellGeom::methodname(const int i) const
{
  switch (i)
  {
  case PHG4CylinderCellDefs::sizebinning:
    return "Bins in cm";
    break;
  case PHG4CylinderCellDefs::etaphibinning:
    return "Eta/Phi bins";
    break;
  case PHG4CylinderCellDefs::etaslatbinning:
    return "Eta/numslat bins";
    break;
  default:
    break;
  }
  return "Unknown";
}

void PHG4BlockCellGeom::check_binning_method(const int i) const
{
  if (_binning != i)
  {
    cout << "different binning method used " << methodname(_binning)
         << ", not : " << methodname(i)
         << endl;
    exit(1);
  }
  return;
}

void PHG4BlockCellGeom::check_binning_method_eta(const std::string& src) const
{
  if (_binning != PHG4CylinderCellDefs::etaphibinning &&
      _binning != PHG4CylinderCellDefs::etaslatbinning)
  {
    if (src.size())
      cout << src << " : ";

    cout << "different binning method used " << methodname(_binning)
         << ", not : " << methodname(PHG4CylinderCellDefs::etaphibinning)
         << " or " << methodname(PHG4CylinderCellDefs::etaslatbinning)
         << endl;
    exit(1);
  }
  return;
}

void PHG4BlockCellGeom::check_binning_method_x(const std::string& src) const
{
  if (_binning != PHG4CylinderCellDefs::etaphibinning &&
      _binning != PHG4CylinderCellDefs::sizebinning &&
      _binning != PHG4CylinderCellDefs::etaslatbinning)
  {
    if (src.size())
      cout << src << " : ";

    cout << "different binning method used " << methodname(_binning)
         << ", not : " << methodname(PHG4CylinderCellDefs::etaphibinning)
         << " or " << methodname(PHG4CylinderCellDefs::sizebinning)
         << endl;
    exit(1);
  }
  return;
}
