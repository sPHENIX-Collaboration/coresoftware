#include "PHField3DCartesian.h"

#include <phool/phool.h>

#include <TDirectory.h>  // for TDirectory, gDirectory
#include <TFile.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <TNtuple.h>
#include <TSystem.h>
#pragma GCC diagnostic pop

#include <Geant4/G4SystemOfUnits.hh>

// stacktrace gives a shadow warning
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/stacktrace.hpp>
#pragma GCC diagnostic pop

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <set>
#include <utility>


PHField3DCartesian::PHField3DCartesian(const std::string &fname, const float magfield_rescale, const float innerradius, const float outerradius, const float size_z)
  : filename(fname)
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
  std::cout << "\n================ Begin Construct Mag Field =====================" << std::endl;
  std::cout << "\n-----------------------------------------------------------"
            << "\n      Magnetic field Module - Verbosity:"
            << "\n-----------------------------------------------------------";

  // open file
  TFile *rootinput = TFile::Open(filename.c_str());
  if (!rootinput)
  {
    std::cout << "\n could not open " << filename << " exiting now" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  std::cout << "\n ---> "
               "Reading the field grid from "
            << filename << " ... " << std::endl;

  //  get root NTuple objects
  TNtuple *field_map = nullptr;
  rootinput->GetObject("fieldmap", field_map);
  if (field_map == nullptr)
  {
    std::cout << PHWHERE << " Could not load fieldmap ntuple from "
	      << filename << " exiting now" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  Float_t ROOT_X, ROOT_Y, ROOT_Z;
  Float_t ROOT_BX, ROOT_BY, ROOT_BZ;
  field_map->SetBranchAddress("x", &ROOT_X);
  field_map->SetBranchAddress("y", &ROOT_Y);
  field_map->SetBranchAddress("z", &ROOT_Z);
  field_map->SetBranchAddress("bx", &ROOT_BX);
  field_map->SetBranchAddress("by", &ROOT_BY);
  field_map->SetBranchAddress("bz", &ROOT_BZ);
  for (int i = 0; i < field_map->GetEntries(); i++)
  {
    field_map->GetEntry(i);
    trio coord_key(ROOT_X * cm, ROOT_Y * cm, ROOT_Z * cm);
    trio field_val(ROOT_BX * tesla * magfield_rescale, ROOT_BY * tesla * magfield_rescale, ROOT_BZ * tesla * magfield_rescale);
    xvals.insert(ROOT_X * cm);
    yvals.insert(ROOT_Y * cm);
    zvals.insert(ROOT_Z * cm);
    if ((std::sqrt(ROOT_X * cm * ROOT_X * cm + ROOT_Y * cm * ROOT_Y * cm) >= innerradius &&
	 std::sqrt(ROOT_X * cm * ROOT_X * cm + ROOT_Y * cm * ROOT_Y * cm) <= outerradius) ||
	std::abs(ROOT_Z * cm) > size_z )
    {
      fieldmap[coord_key] = field_val;
    }
  }
  xmin = *(xvals.begin());
  xmax = *(xvals.rbegin());

  ymin = *(yvals.begin());
  ymax = *(yvals.rbegin());
  if (ymin != xmin || ymax != xmax)
  {
    std::cout << "PHField3DCartesian: Compiler bug!!!!!!!! Do not use inlining!!!!!!" << std::endl;
    std::cout << "exiting now - recompile with -fno-inline" << std::endl;
    exit(1);
  }

  zmin = *(zvals.begin());
  zmax = *(zvals.rbegin());

  xstepsize = (xmax - xmin) / (xvals.size() - 1);
  ystepsize = (ymax - ymin) / (yvals.size() - 1);
  zstepsize = (zmax - zmin) / (zvals.size() - 1);

  delete field_map;
  delete rootinput;
  std::cout << "\n================= End Construct Mag Field ======================\n"
            << std::endl;
}

PHField3DCartesian::~PHField3DCartesian()
{
  if (Verbosity() > 0)
  {
    std::cout << "PHField3DCartesian: cache hits: " << cache_hits
              << " cache misses: " << cache_misses
              << std::endl;
  }
}

void PHField3DCartesian::GetFieldValue(const double point[4], double *Bfield) const
{
  static double xsav = -1000000.;
  static double ysav = -1000000.;
  static double zsav = -1000000.;

  double x = point[0];
  double y = point[1];
  double z = point[2];

  Bfield[0] = 0.0;
  Bfield[1] = 0.0;
  Bfield[2] = 0.0;
  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
  {
    static int ifirst = 0;
    if (ifirst < 10)
    {
      std::cout << "PHField3DCartesian::GetFieldValue: "
                << "Invalid coordinates: "
                << "x: " << x / cm
                << ", y: " << y / cm
                << ", z: " << z / cm
                << " bailing out returning zero bfield"
                << std::endl;
      std::cout << "previous point: "
                << "x: " << xsav / cm
                << ", y: " << ysav / cm
                << ", z: " << zsav / cm
                << std::endl;
      std::cout << "Here is the stacktrace: " << std::endl;
      std::cout << boost::stacktrace::stacktrace();
      std::cout << "This is not a segfault. Check the stacktrace for the guilty party (typically #2)" << std::endl;
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
  double xkey[2];
  std::set<float>::const_iterator it = xvals.lower_bound(x);
  xkey[0] = *it;
  if (it == xvals.begin())
  {
    xkey[1] = *it;
    if (x < xkey[0])
    {
      std::cout << PHWHERE << ": This should not happen! x too small - outside range: " << x / cm << std::endl;
      return;
    }
  }
  else
  {
    --it;
    xkey[1] = *it;
  }

  double ykey[2];
  it = yvals.lower_bound(y);
  ykey[0] = *it;
  if (it == yvals.begin())
  {
    ykey[1] = *it;
    if (y < ykey[0])
    {
      std::cout << PHWHERE << ": This should not happen! y too small - outside range: " << y / cm << std::endl;
      return;
    }
  }
  else
  {
    --it;
    ykey[1] = *it;
  }
  double zkey[2];
  it = zvals.lower_bound(z);
  zkey[0] = *it;
  if (it == zvals.begin())
  {
    zkey[1] = *it;
    if (z < zkey[0])
    {
      std::cout << PHWHERE << ": This should not happen! z too small - outside range: " << z / cm << std::endl;
      return;
    }
  }
  else
  {
    --it;
    zkey[1] = *it;
  }
  if (xkey_save != xkey[0] ||
      ykey_save != ykey[0] ||
      zkey_save != zkey[0])
  {
    cache_misses++;
    xkey_save = xkey[0];
    ykey_save = ykey[0];
    zkey_save = zkey[0];

    std::map<boost::tuple<float, float, float>, boost::tuple<float, float, float> >::const_iterator magval;
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
            std::cout << PHWHERE << " could not locate key in " << filename
		      << " value: x: " << xkey[i] / cm
                      << ", y: " << ykey[j] / cm
                      << ", z: " << zkey[k] / cm << std::endl;
            return;
          }
          xyz[i][j][k][0] = (magval->first).get<0>();
          xyz[i][j][k][1] = (magval->first).get<1>();
          xyz[i][j][k][2] = (magval->first).get<2>();
          bf[i][j][k][0] = (magval->second).get<0>();
          bf[i][j][k][1] = (magval->second).get<1>();
          bf[i][j][k][2] = (magval->second).get<2>();
          if (Verbosity() > 0)
          {
            std::cout << "read x/y/z: " << xyz[i][j][k][0] / cm << "/"
                      << xyz[i][j][k][1] / cm << "/"
                      << xyz[i][j][k][2] / cm << " bx/by/bz: "
                      << bf[i][j][k][0] / tesla << "/"
                      << bf[i][j][k][1] / tesla << "/"
                      << bf[i][j][k][2] / tesla << std::endl;
          }
        }
      }
    }
  }
  else
  {
    cache_hits++;
  }

  // how far are we away from the reference point
  double xinblock = point[0] - xkey[1];
  double yinblock = point[1] - ykey[1];
  double zinblock = point[2] - zkey[1];
  // normalize distance to step size
  double fractionx = xinblock / xstepsize;
  double fractiony = yinblock / ystepsize;
  double fractionz = zinblock / zstepsize;
  if (Verbosity() > 0)
  {
    std::cout << "x/y/z stepsize: " << xstepsize / cm << "/" << ystepsize / cm << "/" << zstepsize / cm << std::endl;
    std::cout << "x/y/z inblock: " << xinblock / cm << "/" << yinblock / cm << "/" << zinblock / cm << std::endl;
    std::cout << "x/y/z fraction: " << fractionx << "/" << fractiony << "/" << fractionz << std::endl;
  }

  // linear extrapolation in cube:

  //Vxyz =
  //V000 * x * y * z +
  //V100 * (1 - x) * y * z +
  //V010 * x * (1 - y) * z +
  //V001 * x y * (1 - z) +
  //V101 * (1 - x) * y * (1 - z) +
  //V011 * x * (1 - y) * (1 - z) +
  //V110 * (1 - x) * (1 - y) * z +
  //V111 * (1 - x) * (1 - y) * (1 - z)

  for (int i = 0; i < 3; i++)
  {
    Bfield[i] = bf[0][0][0][i] * fractionx * fractiony * fractionz +
                bf[1][0][0][i] * (1. - fractionx) * fractiony * fractionz +
                bf[0][1][0][i] * fractionx * (1. - fractiony) * fractionz +
                bf[0][0][1][i] * fractionx * fractiony * (1. - fractionz) +
                bf[1][0][1][i] * (1. - fractionx) * fractiony * (1. - fractionz) +
                bf[0][1][1][i] * fractionx * (1. - fractiony) * (1. - fractionz) +
                bf[1][1][0][i] * (1. - fractionx) * (1. - fractiony) * fractionz +
                bf[1][1][1][i] * (1. - fractionx) * (1. - fractiony) * (1. - fractionz);
  }

  return;
}
