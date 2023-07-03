#ifndef BBC_BBCDEFS_H
#define BBC_BBCDEFS_H

#include <TMath.h>
#include <iostream>

namespace BbcDefs
{
  const Double_t index_refract = 1.4585;
  const Double_t v_ckov = 1.0 / index_refract;  // velocity threshold for CKOV
  const Double_t C = 29.9792458;                // cm/ns


}  // namespace BbcDefs

#endif
