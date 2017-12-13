#ifndef TRACKERDEFS_H
#define TRACKERDEFS_H

namespace TrackerDefs
{
typedef unsigned long long keytype;

// CINT does not know the __attribute__((unused))
#ifndef __CINT__
// __attribute__((unused)) prevents warnings about unused variables
// key layout:
//  common upper 16 bits
//   56 - 64  tracker id
//   48 - 56  layer
static unsigned int bitshift_trackerid __attribute__((unused)) = 64 - 8;
static unsigned int bitshift_layer __attribute__((unused)) = bitshift_trackerid - 8;
//  tracker id specific next 16 bits
//
static unsigned int bitshift_ladder __attribute__((unused)) = bitshift_layer - 8;
static unsigned int bitshift_chip __attribute__((unused)) = bitshift_ladder - 8;
//  implementation dependent lower 32 bits
//
static unsigned int bitshift_row __attribute__((unused)) = 16;
static unsigned int bitshift_col __attribute__((unused)) = 0;
#endif


char get_trackerid(const TrackerDefs::keytype key);
char get_layer(const TrackerDefs::keytype key);
long get_index(const TrackerDefs::keytype key);

/// Enumeration for tracker id to easily maintain consistency
enum TRACKERID
{
  mvtx_id = 0,
  intt_id = 1,
  tpc_id = 2
};


/// MVTX binning used for hits & clusters
namespace MVTXBinning
{
keytype genhitkey(const char trackerid, const char layer,
               const char ladder, const char chip,
               const unsigned int row, const unsigned int col);
keytype gencluskey(const char trackerid, const char layer,
               const char ladder, const char chip,
               const long bit32_index);

char get_ladder(const TrackerDefs::keytype key);
char get_chip(const TrackerDefs::keytype key);
unsigned int get_row(const TrackerDefs::keytype key);
unsigned int get_col(const TrackerDefs::keytype key);
};

/// INTT binning used for hits & clusters
namespace INTTBinning
{

};

/// TPC binning used for hits & clusters
namespace TPCBinning
{

};
}

#endif
