#ifndef CALOTOWERDEFS_H
#define CALOTOWERDEFS_H

namespace CaloTowerDefs
{
  enum DetectorSystem
  {
    CEMC = 0,
    HCALIN = 1,
    HCALOUT = 2,
    SEPD = 3,
    ZDC = 4
  };

  enum BuilderType
  {
    kPRDFTowerv1 = 0,
    kPRDFWaveform = 1,
    kWaveformTowerv2 = 2
  };
}

#endif
