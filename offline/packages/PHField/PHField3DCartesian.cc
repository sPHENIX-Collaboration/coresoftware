#include "PHField3DCartesian.h"

#include <TDirectory.h>                        // for TDirectory, gDirectory
#include <TFile.h>
#include <TNtuple.h>

#include <Geant4/G4SystemOfUnits.hh>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <utility>

using namespace std;

typedef boost::tuple<double, double, double> trio;
std::map<boost::tuple<double, double, double>, boost::tuple<double, double, double> > fieldmap;
std::set<double> xvals;
std::set<double> yvals;
std::set<double> zvals;

PHField3DCartesian::PHField3DCartesian(const string &fname, const float magfield_rescale)
  : filename(fname)
  , xmin(1000000)
  , xmax(-1000000)
  , ymin(1000000)
  , ymax(-1000000)
  , zmin(1000000)
  , zmax(-1000000)
  , xstepsize(NAN)
  , ystepsize(NAN)
  , zstepsize(NAN)
  , xkey_save(NAN)
  , ykey_save(NAN)
  , zkey_save(NAN)
  , cache_hits(0)
  , cache_misses(0)
{
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      for (int k = 0; k < 2; k++)
      {
        for (int l = 0; l < 3; l++)
        {
          xyz[i][j][k][l] = NAN;
          bf[i][j][k][l] = NAN;
        }
      }
    }
  }
  cout << "\n================ Begin Construct Mag Field =====================" << endl;
  cout << "\n-----------------------------------------------------------"
       << "\n      Magnetic field Module - Verbosity:"
       << "\n-----------------------------------------------------------";

  // open file
  TFile *rootinput = TFile::Open(filename.c_str());
  if (!rootinput)
  {
    cout << "\n could not open " << filename << " exiting now" << endl;
    exit(1);
  }
  cout << "\n ---> "
          "Reading the field grid from "
       << filename << " ... " << endl;
  rootinput->cd();

  //  get root NTuple objects
  TNtuple *field_map = (TNtuple *) gDirectory->Get("ntuple");
  Float_t ROOT_X, ROOT_Y, ROOT_Z;
  Float_t ROOT_BX, ROOT_BY, ROOT_BZ;
  field_map->SetBranchAddress("x", &ROOT_X);
  field_map->SetBranchAddress("y", &ROOT_Y);
  field_map->SetBranchAddress("z", &ROOT_Z);
  field_map->SetBranchAddress("Bx", &ROOT_BX);
  field_map->SetBranchAddress("By", &ROOT_BY);
  field_map->SetBranchAddress("Bz", &ROOT_BZ);

  for (int i = 0; i < field_map->GetEntries(); i++)
  {
    field_map->GetEntry(i);
    trio coord_key(ROOT_X * m, ROOT_Y * m, ROOT_Z * m);
    trio field_val(ROOT_BX * tesla, ROOT_BY * tesla, ROOT_BZ * tesla);
    xvals.insert(ROOT_X * m);
    yvals.insert(ROOT_Y * m);
    zvals.insert(ROOT_Z * m);
    fieldmap[coord_key] = field_val;
  }
  xmin = *(xvals.begin());
  xmax = *(xvals.rbegin());

  ymin = *(yvals.begin());
  ymax = *(yvals.rbegin());
  if (ymin != xmin || ymax != xmax)
  {
    cout << "PHField3DCartesian: Compiler bug!!!!!!!! Do not use inlining!!!!!!" << endl;
    cout << "exiting now - recompile with -fno-inline" << endl;
    exit(1);
  }

  if (magfield_rescale != 1.0)
  {
    cout << "PHField3DCartesian: Rescale not implemented" << endl;
    exit(1);
  }

  zmin = *(zvals.begin());
  zmax = *(zvals.rbegin());

  xstepsize = (xmax - xmin) / (xvals.size() - 1);
  ystepsize = (ymax - ymin) / (yvals.size() - 1);
  zstepsize = (zmax - zmin) / (zvals.size() - 1);

  if (rootinput) rootinput->Close();

  cout << "\n================= End Construct Mag Field ======================\n"
       << endl;
}

PHField3DCartesian::~PHField3DCartesian()
{
  //   cout << "PHField3DCartesian: cache hits: " << cache_hits
  //        << " cache misses: " << cache_misses
  //        << endl;
}

void PHField3DCartesian::GetFieldValue(const double point[4], double *Bfield) const
{
  cout << "untested code - I don't know if this is being used, drop me a line (with the field) and I test this"
       << endl
       << "   Chris P." << endl;
  assert(0);

  static double xsav = -1000000.;
  static double ysav = -1000000.;
  static double zsav = -1000000.;

  double x = point[0];
  double y = point[1];
  double z = point[2];

  Bfield[0] = 0.0;
  Bfield[1] = 0.0;
  Bfield[2] = 0.0;
  if (!isfinite(x) || !isfinite(y) || !isfinite(z))
  {
    static int ifirst = 0;
    if (ifirst < 10)
    {
      cout << "PHField3DCartesian::GetFieldValue: "
           << "Invalid coordinates: "
           << "x: " << x
           << ", y: " << y
           << ", z: " << z
           << " bailing out returning zero bfield"
           << endl;
      cout << "previous point: "
           << "x: " << xsav
           << ", y: " << ysav
           << ", z: " << zsav
           << endl;
      ifirst++;
    }
    return;
  }
  xsav = x;
  ysav = y;
  zsav = z;

  if (point[0] < xmin || point[0] > xmax ||
      point[1] < ymin || point[1] > ymax ||
      point[2] < zmin || point[2] > zmax)
  {
    return;
  }
  set<double>::const_iterator it = xvals.lower_bound(x);
  if (it == xvals.begin())
  {
    cout << "x too small - outside range: " << x << endl;
    return;
  }
  double xkey[2];
  xkey[0] = *it;
  --it;
  xkey[1] = *it;

  it = yvals.lower_bound(y);
  if (it == yvals.begin())
  {
    cout << "y too small - outside range: " << y << endl;
    return;
  }
  double ykey[2];
  ykey[0] = *it;
  --it;
  ykey[1] = *it;

  it = zvals.lower_bound(z);
  if (it == zvals.begin())
  {
    cout << "z too small - outside range: " << z << endl;
    return;
  }
  double zkey[2];
  zkey[0] = *it;
  --it;
  zkey[1] = *it;

  if (xkey_save != xkey[0] ||
      ykey_save != ykey[0] ||
      zkey_save != zkey[0])
  {
    cache_misses++;
    xkey_save = xkey[0];
    ykey_save = ykey[0];
    zkey_save = zkey[0];

    map<boost::tuple<double, double, double>, boost::tuple<double, double, double> >::const_iterator magval;
    trio key;
    for (int i = 0; i < 2; i++)
    {
      for (int j = 0; j < 2; j++)
      {
        for (int k = 0; k < 2; k++)
        {
          key = boost::make_tuple(xkey[i], ykey[j], zkey[k]);
          magval = fieldmap.find(key);
          if (magval == fieldmap.end())
          {
            cout << "could not locate key, x: " << xkey[i] / m
                 << ", y: " << ykey[j] / m
                 << ", z: " << zkey[k] / m << endl;
            return;
          }
          xyz[i][j][k][0] = (magval->first).get<0>();
          xyz[i][j][k][1] = (magval->first).get<1>();
          xyz[i][j][k][2] = (magval->first).get<2>();
          bf[i][j][k][0] = (magval->second).get<0>();
          bf[i][j][k][1] = (magval->second).get<1>();
          bf[i][j][k][2] = (magval->second).get<2>();
          cout << "read x/y/z: " << xyz[i][j][k][0] << "/"
               << xyz[i][j][k][1] << "/"
               << xyz[i][j][k][2] << " bx/by/bz: "
               << bf[i][j][k][0] << "/"
               << bf[i][j][k][1] << "/"
               << bf[i][j][k][2] << endl;
        }
      }
    }
  }
  else
  {
    cache_hits++;
  }
  // linear extrapolation in cube:
  //Vxyz = 	V000 (1 - x) (1 - y) (1 - z) +
  //V100 x (1 - y) (1 - z) +
  //V010 (1 - x) y (1 - z) +
  //V001 (1 - x) (1 - y) z +
  //V101 x (1 - y) z +
  //V011 (1 - x) y z +
  //V110 x y (1 - z) +
  //V111 x y z

  double xinblock = point[0] - xkey[1];
  double yinblock = point[1] - ykey[1];
  double zinblock = point[2] - zkey[1];
  //   cout << "x/y/z stepsize: " << xstepsize << "/" << ystepsize << "/" << zstepsize << endl;
  //   cout << "x/y/z inblock: " << xinblock << "/" << yinblock << "/" << zinblock << endl;
  for (int i = 0; i < 3; i++)
  {
    Bfield[i] = bf[0][0][0][i] * (xstepsize - xinblock) * (ystepsize - yinblock) * (zstepsize - zinblock) +
                bf[1][0][0][i] * xinblock * (ystepsize - yinblock) * (zstepsize - zinblock) +
                bf[0][1][0][i] * (xstepsize - xinblock) * yinblock * (zstepsize - zinblock) +
                bf[0][0][1][i] * (xstepsize - xinblock) * (ystepsize - yinblock) * zinblock +
                bf[1][0][1][i] * xinblock * (ystepsize - yinblock) * zinblock +
                bf[0][1][1][i] * (xstepsize - xinblock) * yinblock * zinblock +
                bf[1][1][0][i] * xinblock * yinblock * (zstepsize - zinblock) +
                bf[1][1][1][i] * xinblock * yinblock * zinblock;
  }

  return;
}
