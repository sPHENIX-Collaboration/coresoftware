#ifndef __TrackerDefs_H__
#define __TrackerDefs_H__

#ifdef __CINT__
#include <stdint.h>
#else
#include <cstdint>
#endif
#include <limits.h>
// #include <climits>
#include <iostream>

namespace TrackerDefs
{
typedef uint32_t hitkeytype; // 32 bit TrackerHit id type
typedef uint64_t keytype;    // 64 but TrackerCluster id type

// CINT does not know the __attribute__((unused))
// __attribute__((unused)) prevents warnings about unused variables
#ifndef __CINT__

// D. McGlinchey - I don't understand why this needs to be hidden from
//                 the dictionary generation ...
static hitkeytype HITKEYMAX = ULONG_MAX;
static keytype KEYMAX = ULLONG_MAX;

// key layout:
//  common upper 16 bits
//   56 - 64  tracker id
//   48 - 56  layer
static unsigned int bitshift_trackerid __attribute__((unused)) = 32 - 8;
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

/// Get the tracker ID from either key type
char get_trackerid(const TrackerDefs::hitkeytype key);
char get_trackerid(const TrackerDefs::keytype key);

/// Get the layer number from either key type
char get_layer(const TrackerDefs::hitkeytype key);
char get_layer(const TrackerDefs::keytype key);

/// Get the lower 32 bits for cluster keys only
long get_index(const TrackerDefs::keytype key);

/// Print the bits for either key type
void print_bits(const TrackerDefs::hitkeytype key, std::ostream& os = std::cout);
void print_bits(const TrackerDefs::keytype key, std::ostream& os = std::cout);

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
hitkeytype genhitkey(const char trackerid, const char layer,
               const char ladder, const char chip);
keytype gencluskey(const char trackerid, const char layer,
               const char ladder, const char chip,
               const unsigned int bit32_index);

char get_ladder(const TrackerDefs::hitkeytype key);
char get_ladder(const TrackerDefs::keytype key);

char get_chip(const TrackerDefs::hitkeytype key);
char get_chip(const TrackerDefs::keytype key);

unsigned short get_row(const TrackerDefs::keytype key);
unsigned short get_col(const TrackerDefs::keytype key);
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
