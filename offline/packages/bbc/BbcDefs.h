#ifndef __BBC_BBCDEFS_H__
#define __BBC_BBCDEFS_H__

#include <gsl/gsl_const_cgsm.h>

namespace BbcDefs
{
  const double index_refract = 1.4585;        // quartz radiator index of refraction
  const double v_ckov = 1.0 / index_refract;  // velocity threshold for CKOV
  //const double C = 29.9792458;                // cm/ns
  const double C = GSL_CONST_CGSM_SPEED_OF_LIGHT*1e-9;                // cm/ns

  const int BBC_N_PMT = 128;
  const int BBC_N_FEECH = 256;
  const int MAX_SAMPLES = 31;

}  // namespace BbcDefs

#endif  // __BBC_BBCDEFS_H__
