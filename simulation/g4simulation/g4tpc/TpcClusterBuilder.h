#ifndef G4TPC_TPCCLUSTERBUILDER_H
#define G4TPC_TPCCLUSTERBUILDER_H

#include <trackbase/ActsGeometry.h>
#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>
#include <map>
#include <climits>

class TrkrCluster;
class PHG4TpcCylinderGeom;
  // This structure collects data from PHG$TpcPadPlaneReadout for cluster data in a given layer for a given truth particle.

class TpcClusterBuilder 
{

  using PairCluskeyCluster = std::pair<TrkrDefs::cluskey,TrkrCluster*>;

  public:
  static constexpr double AdcClockPeriod = 53.0; // ns (copied from TpcClusterizer.h)
  static ActsGeometry* tGeometry;

  TrkrDefs::hitsetkey   hitsetkey { UCHAR_MAX };
  PHG4TpcCylinderGeom *layerGeom { nullptr }; // unique to the layer
  bool   has_data       { false };
  short  layer          {  SHRT_MIN };
  unsigned int side     { UINT_MAX };
  int    neff_electrons { 0  };
  double phi_integral   { 0. };
  double time_integral  { 0. };
  int    phi_bin_lo     { INT_MAX };
  int    phi_bin_hi     { INT_MIN };
  int    time_bin_lo    { INT_MAX };
  int    time_bin_hi    { INT_MIN };
  int    nphibins       { INT_MIN };
  bool   hasPhiBins     { false   };
  bool   hasTimeBins    { false   };

  TpcClusterBuilder& operator+=(const TpcClusterBuilder& rhs);

  void fillPhiBins  (const std::vector<int>& bins);
  void fillTimeBins (const std::vector<int>& bins);
  void set_has_data();

  TrkrCluster* build() const;

  TpcClusterBuilder() {};
  ~TpcClusterBuilder(){};

};

#endif  //TRACKBASE_PADPLANEREADOUTSTRUCT_H
