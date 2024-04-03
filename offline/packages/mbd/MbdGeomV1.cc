#include "MbdGeomV1.h"

#include <cmath>

// kludge where we have the hardcoded positions of the tubes
// Should really be put in a database
// These are the x,y for the south MBD (in cm).
// The north inverts the x coordinate (x -> -x)
static const float PmtLoc[64][2] = {
    {-12.2976, 4.26},
    {-12.2976, 1.42},
    {-9.83805, 8.52},
    {-9.83805, 5.68},
    {-9.83805, 2.84},
    {-7.37854, 9.94},
    {-7.37854, 7.1},
    {-7.37854, 4.26},
    {-7.37854, 1.42},
    {-4.91902, 11.36},
    {-4.91902, 8.52},
    {-4.91902, 5.68},
    {-2.45951, 12.78},
    {-2.45951, 9.94},
    {-2.45951, 7.1},
    {0, 11.36},
    {0, 8.52},
    {2.45951, 12.78},
    {2.45951, 9.94},
    {2.45951, 7.1},
    {4.91902, 11.36},
    {4.91902, 8.52},
    {4.91902, 5.68},
    {7.37854, 9.94},
    {7.37854, 7.1},
    {7.37854, 4.26},
    {7.37854, 1.42},
    {9.83805, 8.52},
    {9.83805, 5.68},
    {9.83805, 2.84},
    {12.2976, 4.26},
    {12.2976, 1.42},
    {12.2976, -4.26},
    {12.2976, -1.42},
    {9.83805, -8.52},
    {9.83805, -5.68},
    {9.83805, -2.84},
    {7.37854, -9.94},
    {7.37854, -7.1},
    {7.37854, -4.26},
    {7.37854, -1.42},
    {4.91902, -11.36},
    {4.91902, -8.52},
    {4.91902, -5.68},
    {2.45951, -12.78},
    {2.45951, -9.94},
    {2.45951, -7.1},
    {0, -11.36},
    {0, -8.52},
    {-2.45951, -12.78},
    {-2.45951, -9.94},
    {-2.45951, -7.1},
    {-4.91902, -11.36},
    {-4.91902, -8.52},
    {-4.91902, -5.68},
    {-7.37854, -9.94},
    {-7.37854, -7.1},
    {-7.37854, -4.26},
    {-7.37854, -1.42},
    {-9.83805, -8.52},
    {-9.83805, -5.68},
    {-9.83805, -2.84},
    {-12.2976, -4.26},
    {-12.2976, -1.42}};

MbdGeomV1::MbdGeomV1()
{
  // Set the pmt locations
  for (unsigned int ipmt = 0; ipmt < 128; ipmt++)
  {
    int arm = ipmt / 64;

    float xsign = 1.;
    float zsign = -1.;
    if (arm == 1)  // north
    {
      xsign = -1.;
      zsign = 1.;
    }

    float tube_x = xsign * PmtLoc[ipmt % 64][0];
    float tube_y = PmtLoc[ipmt % 64][1];
    float tube_z = zsign * 253.;

    MbdGeomV1::set_xyz(ipmt, tube_x, tube_y, tube_z);
  }
}

void MbdGeomV1::set_xyz(const unsigned int ipmt, const float x, const float y, const float z)
{
  pmt_x[ipmt] = x;
  pmt_y[ipmt] = y;
  pmt_z[ipmt] = z;
  pmt_r[ipmt] = std::sqrt(x * x + y * y);
  pmt_phi[ipmt] = std::atan2(y, x);
}

