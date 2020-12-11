#include "PHFieldCleo.h"

#include <Geant4/G4SystemOfUnits.hh>

#include <TSystem.h>

#include <cmath>                        // for modf, fabs, isfinite, NAN
#include <cstdlib>                     // for exit
#include <fstream>
#include <iostream>
#include <memory>                       // for allocator_traits<>::value_type

using namespace std;

PHFieldCleo::PHFieldCleo(const string &filename, const int verb, const float magfield_rescale)
  : PHField(verb)
  , m_MagFieldScale(magfield_rescale)
{
  ifstream file(filename);
  if (!file.is_open())
  {
    cout << "Could not open " << filename << " exiting now" << endl;
    gSystem->Exit(1);
    exit(1);
  }
  char buffer[256];
  // Read table dimensions
  int tmp;
  file >> nz >> ny >> nx >> tmp;  // first z, then y and x (y,x are 50, it does not matter)
  xField.resize(nx);
  yField.resize(nx);
  zField.resize(nx);
  int ix, iy, iz;
  for (ix = 0; ix < nx; ix++)
  {
    xField[ix].resize(ny);
    yField[ix].resize(ny);
    zField[ix].resize(ny);
    for (iy = 0; iy < ny; iy++)
    {
      xField[ix][iy].resize(nz);
      yField[ix][iy].resize(nz);
      zField[ix][iy].resize(nz);
    }
  }

  // Ignore other header information
  // The first line whose second character is '0' is considered to
  // be the last line of the header.
  do
  {
    file.getline(buffer, 256);
  } while (buffer[1] != '0');

  // Read in the data
  double xval, yval, zval, bx, by, bz;
  double save_x = NAN;
  double save_y = NAN;
  for (ix = 0; ix < nx; ix++)
  {
    for (iy = 0; iy < ny; iy++)
    {
      for (iz = 0; iz < nz; iz++)
      {
        file >> xval >> yval >> zval >> bx >> by >> bz;
        if (!isfinite(save_x))
        {
          save_x = xval;
        }
        else
        {
          if (save_x != xval)
          {
            cout << "mismatch expected x coordinate: " << save_x
                 << " read " << xval << endl;
            gSystem->Exit(1);
          }
        }
        if (!isfinite(save_y))
        {
          save_y = yval;
        }
        else
        {
          if (save_y != yval)
          {
            cout << "mismatch expected y coordinate: " << save_y
                 << " read " << yval << endl;
            gSystem->Exit(1);
          }
        }
        if (ix == 0 && iy == 0 && iz == 0)
        {
          minx = xval * cm;
          miny = yval * cm;
          minz = zval * cm;
        }
        xField[ix][iy][iz] = bx * gauss * m_MagFieldScale;
        yField[ix][iy][iz] = by * gauss * m_MagFieldScale;
        zField[ix][iy][iz] = bz * gauss * m_MagFieldScale;
      }
      save_y = NAN;
    }
    save_x = NAN;
  }
  file.close();
  maxx = xval * cm;
  maxy = yval * cm;
  maxz = zval * cm;
  dx = maxx - minx;
  dy = maxy - miny;
  dz = maxz - minz;
}

void PHFieldCleo::GetFieldValue(const double point[4], double *Bfield) const
{
  double x = fabs(point[0]);
  double y = fabs(point[1]);
  double z = point[2];

  // Check that the point is within the defined region
  if (x >= minx && x <= maxx && y >= miny && y <= maxy && z >= minz && z <= maxz)
  {
    // Position of given point within region, normalized to the range
    // [0,1]
    double xfraction = (x - minx) / dx;
    double yfraction = (y - miny) / dy;
    double zfraction = (z - minz) / dz;

    // Need addresses of these to pass to modf below.
    // modf uses its second argument as an OUTPUT argument.
    double xdindex, ydindex, zdindex;

    // Position of the point within the cuboid defined by the
    // nearest surrounding tabulated points
    double xlocal = (std::modf(xfraction * (nx - 1), &xdindex));
    double ylocal = (std::modf(yfraction * (ny - 1), &ydindex));
    double zlocal = (std::modf(zfraction * (nz - 1), &zdindex));

    // The indices of the nearest tabulated point whose coordinates
    // are all less than those of the given point
    int xindex = static_cast<int>(xdindex);
    int yindex = static_cast<int>(ydindex);
    int zindex = static_cast<int>(zdindex);

    // Full 3-dimensional version
    Bfield[0] = xField[xindex][yindex][zindex] * (1 - xlocal) * (1 - ylocal) * (1 - zlocal) + xField[xindex][yindex][zindex + 1] * (1 - xlocal) * (1 - ylocal) * zlocal + xField[xindex][yindex + 1][zindex] * (1 - xlocal) * ylocal * (1 - zlocal) + xField[xindex][yindex + 1][zindex + 1] * (1 - xlocal) * ylocal * zlocal + xField[xindex + 1][yindex][zindex] * xlocal * (1 - ylocal) * (1 - zlocal) + xField[xindex + 1][yindex][zindex + 1] * xlocal * (1 - ylocal) * zlocal + xField[xindex + 1][yindex + 1][zindex] * xlocal * ylocal * (1 - zlocal) + xField[xindex + 1][yindex + 1][zindex + 1] * xlocal * ylocal * zlocal;
    Bfield[1] = yField[xindex][yindex][zindex] * (1 - xlocal) * (1 - ylocal) * (1 - zlocal) + yField[xindex][yindex][zindex + 1] * (1 - xlocal) * (1 - ylocal) * zlocal + yField[xindex][yindex + 1][zindex] * (1 - xlocal) * ylocal * (1 - zlocal) + yField[xindex][yindex + 1][zindex + 1] * (1 - xlocal) * ylocal * zlocal + yField[xindex + 1][yindex][zindex] * xlocal * (1 - ylocal) * (1 - zlocal) + yField[xindex + 1][yindex][zindex + 1] * xlocal * (1 - ylocal) * zlocal + yField[xindex + 1][yindex + 1][zindex] * xlocal * ylocal * (1 - zlocal) + yField[xindex + 1][yindex + 1][zindex + 1] * xlocal * ylocal * zlocal;
    Bfield[2] = zField[xindex][yindex][zindex] * (1 - xlocal) * (1 - ylocal) * (1 - zlocal) + zField[xindex][yindex][zindex + 1] * (1 - xlocal) * (1 - ylocal) * zlocal + zField[xindex][yindex + 1][zindex] * (1 - xlocal) * ylocal * (1 - zlocal) + zField[xindex][yindex + 1][zindex + 1] * (1 - xlocal) * ylocal * zlocal + zField[xindex + 1][yindex][zindex] * xlocal * (1 - ylocal) * (1 - zlocal) + zField[xindex + 1][yindex][zindex + 1] * xlocal * (1 - ylocal) * zlocal + zField[xindex + 1][yindex + 1][zindex] * xlocal * ylocal * (1 - zlocal) + zField[xindex + 1][yindex + 1][zindex + 1] * xlocal * ylocal * zlocal;
  }
  else
  {
    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;
  }
}
