#ifndef G4TPC_TPCCLUSTERBUILDER_H
#define G4TPC_TPCCLUSTERBUILDER_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>
#include <map>
#include <climits>

class TrkrCluster;
class PHG4TpcCylinderGeom;
  // This is really just a structure used in simulation/g4simulation/g4tpc/PHG4TpcElectronDrift.cc
  // to collect necessary statistics from simulation/g4simulation/g4tpc/PHG4TpcPadPlaneReadout.cc
  // 
  // It could have been just std::pair< std::array<int,6>, std::pair<double,double>>, but 
  // having the structure to name those subtypes.


class TpcClusterBuilder 
{
  using MapHitsetkeyUInt   = std::map<TrkrDefs::hitsetkey, unsigned int>;
  using PairCluskeyCluster = std::pair<TrkrDefs::cluskey,TrkrCluster*>;
  public:
    PHG4TpcCylinderGeom *layerGeom { nullptr }; // unique to the layer
    bool   has_data       { false };
    short  layer          {  SHRT_MAX };
    unsigned int side     { 0       };
    int    neff_electrons { 0       };
    double phi_integral   { 0.      };
    double time_integral  { 0.      };
    int    phi_bin_lo     { INT_MAX };
    int    phi_bin_hi     { INT_MIN };
    int    time_bin_lo    { INT_MAX };
    int    time_bin_hi    { INT_MIN };
    int    nphibins       { 0       };
    bool   hasPhiBins     { false   };
    bool   hasTimeBins    { false   };

    TpcClusterBuilder( short _layer,  unsigned int _side, 
        int _neff_electrons, double _phi_integral, 
        double _time_integral, int _phi_bin_lo,  int _phi_bin_hi, 
        int _time_bin_lo, int _time_bin_hi) 
      : layerGeom      { nullptr         }
      , layer          { _layer          }
      , side           { _side           }
      , neff_electrons { _neff_electrons }
      , phi_integral   { _phi_integral   }
      , time_integral  { _time_integral  }
      , phi_bin_lo     { _phi_bin_lo     }
      , phi_bin_hi     { _phi_bin_hi     }
      , time_bin_lo    { _time_bin_lo    }
      , time_bin_hi    { _time_bin_hi    }
    {};
    TpcClusterBuilder() {};
      /* layerGeom {nullptr}, layer {0}, side {0}, neff_electrons {0}, */ 
      /* phi_integral {0.}, time_integral {0.}, */
      /* phi_bin_lo {0}, phi_bin_hi {0}, time_bin_lo {0}, time_bin_hi {0}, nphibins{0} */
    /* {}; */

    TpcClusterBuilder& operator+=(const TpcClusterBuilder& rhs);
    /* void fillBinsLowHigh (int& bin_lo, int&bin_hi, std::vector<int> bins); */
    void fillPhiBins  (const std::vector<int>& bins);
    void fillTimeBins (const std::vector<int>& bins);
    void reset();
    void set_has_data();
    PairCluskeyCluster build(MapHitsetkeyUInt& cluster_cnt) const;

    ~TpcClusterBuilder(){};
};

#endif  //TRACKBASE_PADPLANEREADOUTSTRUCT_H
