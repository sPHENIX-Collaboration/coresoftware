#include "MbdGeomV2.h"

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

MbdGeomV2::MbdGeomV2()
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

    MbdGeomV2::set_xyz(ipmt, tube_x, tube_y, tube_z);
  }

}

void MbdGeomV2::set_xyz(const unsigned int ipmt, const float x, const float y, const float z)
{
  pmt_x[ipmt] = x;
  pmt_y[ipmt] = y;
  pmt_z[ipmt] = z;
  pmt_r[ipmt] = std::sqrt((x * x) + (y * y));
  pmt_phi[ipmt] = std::atan2(y, x);
}

const std::multimap<int,int>& MbdGeomV2::get_hvmap()
{
  if ( pmt_hv.empty() )
  {
    download_hv();
  }
  return pmt_hv;
}

void MbdGeomV2::download_hv()
{
  // South HV Map
  pmt_hv.insert(std::pair<int, int>( 0, 0 ));
  pmt_hv.insert(std::pair<int, int>( 0, 1 ));
  pmt_hv.insert(std::pair<int, int>( 0, 2 ));
  pmt_hv.insert(std::pair<int, int>( 0, 3 ));
  pmt_hv.insert(std::pair<int, int>( 0, 4 ));
  pmt_hv.insert(std::pair<int, int>( 0, 5 ));
  pmt_hv.insert(std::pair<int, int>( 0, 6 ));
  pmt_hv.insert(std::pair<int, int>( 0, 10 ));
  pmt_hv.insert(std::pair<int, int>( 1, 9 ));
  pmt_hv.insert(std::pair<int, int>( 1, 12 ));
  pmt_hv.insert(std::pair<int, int>( 1, 13 ));
  pmt_hv.insert(std::pair<int, int>( 1, 15 ));
  pmt_hv.insert(std::pair<int, int>( 1, 17 ));
  pmt_hv.insert(std::pair<int, int>( 1, 18 ));
  pmt_hv.insert(std::pair<int, int>( 1, 20 ));
  pmt_hv.insert(std::pair<int, int>( 1, 23 ));
  pmt_hv.insert(std::pair<int, int>( 2, 7 ));
  pmt_hv.insert(std::pair<int, int>( 2, 8 ));
  pmt_hv.insert(std::pair<int, int>( 2, 11 ));
  pmt_hv.insert(std::pair<int, int>( 2, 14 ));
  pmt_hv.insert(std::pair<int, int>( 2, 16 ));
  pmt_hv.insert(std::pair<int, int>( 2, 19 ));
  pmt_hv.insert(std::pair<int, int>( 2, 22 ));
  pmt_hv.insert(std::pair<int, int>( 2, 25 ));
  pmt_hv.insert(std::pair<int, int>( 2, 26 ));
  pmt_hv.insert(std::pair<int, int>( 3, 21 ));
  pmt_hv.insert(std::pair<int, int>( 3, 24 ));
  pmt_hv.insert(std::pair<int, int>( 3, 27 ));
  pmt_hv.insert(std::pair<int, int>( 3, 28 ));
  pmt_hv.insert(std::pair<int, int>( 3, 29 ));
  pmt_hv.insert(std::pair<int, int>( 3, 30 ));
  pmt_hv.insert(std::pair<int, int>( 3, 31 ));
  pmt_hv.insert(std::pair<int, int>( 4, 32 ));
  pmt_hv.insert(std::pair<int, int>( 4, 33 ));
  pmt_hv.insert(std::pair<int, int>( 4, 34 ));
  pmt_hv.insert(std::pair<int, int>( 4, 35 ));
  pmt_hv.insert(std::pair<int, int>( 4, 36 ));
  pmt_hv.insert(std::pair<int, int>( 4, 37 ));
  pmt_hv.insert(std::pair<int, int>( 4, 38 ));
  pmt_hv.insert(std::pair<int, int>( 4, 42 ));
  pmt_hv.insert(std::pair<int, int>( 5, 41 ));
  pmt_hv.insert(std::pair<int, int>( 5, 44 ));
  pmt_hv.insert(std::pair<int, int>( 5, 45 ));
  pmt_hv.insert(std::pair<int, int>( 5, 47 ));
  pmt_hv.insert(std::pair<int, int>( 5, 49 ));
  pmt_hv.insert(std::pair<int, int>( 5, 50 ));
  pmt_hv.insert(std::pair<int, int>( 5, 52 ));
  pmt_hv.insert(std::pair<int, int>( 6, 39 ));
  pmt_hv.insert(std::pair<int, int>( 6, 40 ));
  pmt_hv.insert(std::pair<int, int>( 6, 43 ));
  pmt_hv.insert(std::pair<int, int>( 6, 46 ));
  pmt_hv.insert(std::pair<int, int>( 6, 48 ));
  pmt_hv.insert(std::pair<int, int>( 6, 51 ));
  pmt_hv.insert(std::pair<int, int>( 6, 54 ));
  pmt_hv.insert(std::pair<int, int>( 6, 57 ));
  pmt_hv.insert(std::pair<int, int>( 6, 58 ));
  pmt_hv.insert(std::pair<int, int>( 6, 60 ));
  pmt_hv.insert(std::pair<int, int>( 7, 53 ));
  pmt_hv.insert(std::pair<int, int>( 7, 55 ));
  pmt_hv.insert(std::pair<int, int>( 7, 56 ));
  pmt_hv.insert(std::pair<int, int>( 7, 59 ));
  pmt_hv.insert(std::pair<int, int>( 7, 61 ));
  pmt_hv.insert(std::pair<int, int>( 7, 62 ));
  pmt_hv.insert(std::pair<int, int>( 7, 63 ));

  // North HV Map
  pmt_hv.insert(std::pair<int, int>( 8, 64 ));
  pmt_hv.insert(std::pair<int, int>( 8, 65 ));
  pmt_hv.insert(std::pair<int, int>( 8, 66 ));
  pmt_hv.insert(std::pair<int, int>( 8, 67 ));
  pmt_hv.insert(std::pair<int, int>( 8, 68 ));
  pmt_hv.insert(std::pair<int, int>( 8, 70 ));
  pmt_hv.insert(std::pair<int, int>( 8, 74 ));
  pmt_hv.insert(std::pair<int, int>( 8, 77 ));
  pmt_hv.insert(std::pair<int, int>( 9, 76 ));
  pmt_hv.insert(std::pair<int, int>( 9, 81 ));
  pmt_hv.insert(std::pair<int, int>( 9, 84 ));
  pmt_hv.insert(std::pair<int, int>( 9, 87 ));
  pmt_hv.insert(std::pair<int, int>( 9, 91 ));
  pmt_hv.insert(std::pair<int, int>( 9, 94 ));
  pmt_hv.insert(std::pair<int, int>( 9, 95 ));
  pmt_hv.insert(std::pair<int, int>( 10, 71 ));
  pmt_hv.insert(std::pair<int, int>( 10, 72 ));
  pmt_hv.insert(std::pair<int, int>( 10, 75 ));
  pmt_hv.insert(std::pair<int, int>( 10, 78 ));
  pmt_hv.insert(std::pair<int, int>( 10, 80 ));
  pmt_hv.insert(std::pair<int, int>( 10, 83 ));
  pmt_hv.insert(std::pair<int, int>( 10, 86 ));
  pmt_hv.insert(std::pair<int, int>( 10, 89 ));
  pmt_hv.insert(std::pair<int, int>( 10, 90 ));
  pmt_hv.insert(std::pair<int, int>( 11, 69 ));
  pmt_hv.insert(std::pair<int, int>( 11, 73 ));
  pmt_hv.insert(std::pair<int, int>( 11, 79 ));
  pmt_hv.insert(std::pair<int, int>( 11, 82 ));
  pmt_hv.insert(std::pair<int, int>( 11, 85 ));
  pmt_hv.insert(std::pair<int, int>( 11, 88 ));
  pmt_hv.insert(std::pair<int, int>( 11, 92 ));
  pmt_hv.insert(std::pair<int, int>( 11, 93 ));
  pmt_hv.insert(std::pair<int, int>( 12, 96 ));
  pmt_hv.insert(std::pair<int, int>( 12, 97 ));
  pmt_hv.insert(std::pair<int, int>( 12, 98 ));
  pmt_hv.insert(std::pair<int, int>( 12, 99 ));
  pmt_hv.insert(std::pair<int, int>( 12, 100 ));
  pmt_hv.insert(std::pair<int, int>( 12, 102 ));
  pmt_hv.insert(std::pair<int, int>( 12, 106 ));
  pmt_hv.insert(std::pair<int, int>( 12, 108 ));
  pmt_hv.insert(std::pair<int, int>( 12, 109 ));
  pmt_hv.insert(std::pair<int, int>( 13, 101 ));
  pmt_hv.insert(std::pair<int, int>( 13, 105 ));
  pmt_hv.insert(std::pair<int, int>( 13, 111 ));
  pmt_hv.insert(std::pair<int, int>( 13, 119 ));
  pmt_hv.insert(std::pair<int, int>( 14, 103 ));
  pmt_hv.insert(std::pair<int, int>( 14, 104 ));
  pmt_hv.insert(std::pair<int, int>( 14, 107 ));
  pmt_hv.insert(std::pair<int, int>( 14, 110 ));
  pmt_hv.insert(std::pair<int, int>( 14, 112 ));
  pmt_hv.insert(std::pair<int, int>( 14, 115 ));
  pmt_hv.insert(std::pair<int, int>( 14, 116 ));
  pmt_hv.insert(std::pair<int, int>( 14, 118 ));
  pmt_hv.insert(std::pair<int, int>( 14, 121 ));
  pmt_hv.insert(std::pair<int, int>( 14, 122 ));
  pmt_hv.insert(std::pair<int, int>( 15, 113 ));
  pmt_hv.insert(std::pair<int, int>( 15, 114 ));
  pmt_hv.insert(std::pair<int, int>( 15, 117 ));
  pmt_hv.insert(std::pair<int, int>( 15, 120 ));
  pmt_hv.insert(std::pair<int, int>( 15, 123 ));
  pmt_hv.insert(std::pair<int, int>( 15, 124 ));
  pmt_hv.insert(std::pair<int, int>( 15, 125 ));
  pmt_hv.insert(std::pair<int, int>( 15, 126 ));
  pmt_hv.insert(std::pair<int, int>( 15, 127 ));

}

