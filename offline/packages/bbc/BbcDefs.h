#ifndef __BBC_BBCDEFS_H__
#define __BBC_BBCDEFS_H__

#include <TMath.h>
#include <iostream>

namespace BbcDefs
{
  const Double_t index_refract = 1.4585;
  const Double_t v_ckov = 1.0 / index_refract;  // velocity threshold for CKOV
  const Double_t C = 29.9792458;                // cm/ns

  const int BBC_N_PMT = 128;

}  // namespace BbcDefs

#endif  // __BBC_BBCDEFS_H__
