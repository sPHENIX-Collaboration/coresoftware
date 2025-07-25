#ifndef CALOCDB_GEOMETRYCONSTANTS_H
#define CALOCDB_GEOMETRYCONSTANTS_H

namespace CaloGeometry
{
  // EMCAL
  constexpr int CEMC_ETA_BINS = 96;
  constexpr int CEMC_PHI_BINS = 256;
  constexpr int CEMC_NSECTOR = 64;
  constexpr int CEMC_NCHANNEL_PER_SECTOR = 384;
  constexpr int CEMC_NCHANNEL_PER_IB = 64;
  constexpr int CEMC_NIB = 384;
  constexpr int CEMC_NIB_PER_SECTOR = 6;
  constexpr int CEMC_NTOW_IB_SIDE = 8;

  // HCAL
  constexpr int HCAL_ETA_BINS = 24;
  constexpr int HCAL_PHI_BINS = 64;
}  // namespace CaloGeometry

#endif
