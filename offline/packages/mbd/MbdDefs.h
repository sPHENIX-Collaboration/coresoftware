#ifndef __MBD_MBDDEFS_H__
#define __MBD_MBDDEFS_H__

#ifndef ONLINE
#include <gsl/gsl_const_cgsm.h>
#endif

namespace MbdDefs
{
  const double index_refract = 1.4585;        // quartz radiator index of refraction
  const double v_ckov = 1.0 / index_refract;  // velocity threshold for CKOV
#ifndef ONLINE
  const double C = GSL_CONST_CGSM_SPEED_OF_LIGHT*1e-9;                // cm/ns
#else
  const double C = 29.9792458;                // cm/ns
#endif

  const int MBD_N_PMT = 128;
  const int MBD_N_FEECH = 256;
  const int BBC_N_PMT = 128;
  const int BBC_N_FEECH = 256;
  const int MAX_SAMPLES = 31;

}  // namespace MbdDefs

#endif  // __MBD_MBDDEFS_H__
