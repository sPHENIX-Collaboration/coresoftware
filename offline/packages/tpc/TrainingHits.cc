#include "TrainingHits.h"

TrainingHits::TrainingHits()
{
  Reset();
}

void TrainingHits::Reset()
{
  v_adc.fill(0);
  radius = 0.;
  phi = 0.;
  z = 0.;
  phistep = 0.;
  zstep = 0.;
  layer = 0;
  ntouch = 0;
  nedge = 0;
  cluskey = 0;
}
