#ifndef TRACKERDEFS_H
#define TRACKERDEFS_H

namespace TrackerDefs
{
typedef unsigned long long keytype;

// CINT does not know the __attribute__((unused))
#ifndef __CINT__
// __attribute__((unused)) prevents warnings about unused variables
// key layout:
//  common upper 32 bits
//   56 - 64  detector id
//   48 - 56  layer
//   40 - 48  ladder
//   32 - 40  chip
//
static unsigned int bitshift_detid __attribute__((unused)) = 64 - 8;
static unsigned int bitshift_layer __attribute__((unused)) = bitshift_detid - 8;
static unsigned int bitshift_ladder __attribute__((unused)) = bitshift_layer - 8;
static unsigned int bitshift_chip __attribute__((unused)) = bitshift_ladder - 8;
//  implementation dependent lower 32 bits
//
#endif

keytype genkey(const char detid, const char layer,
               const char ladder, const char chip,
               const unsigned int bit32_index);

char get_detid(const TrackerDefs::keytype key);
char get_layer(const TrackerDefs::keytype key);
char get_ladder(const TrackerDefs::keytype key);
char get_chip(const TrackerDefs::keytype key);
unsigned int get_index(const TrackerDefs::keytype key);
}

#endif
