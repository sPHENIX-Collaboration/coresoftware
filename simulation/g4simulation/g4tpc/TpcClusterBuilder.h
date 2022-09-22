#ifndef G4TPC_TPCCLUSTERBUILDER_H
#define G4TPC_TPCCLUSTERBUILDER_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>
#include <map>

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
    PHG4TpcCylinderGeom *layerGeom ; // unique to the layer
    short  layer;
    unsigned int side;
    int    neff_electrons ;
    double phi_integral   ;
    double time_integral  ;
    int    phi_bin_lo     ;
    int    phi_bin_hi     ;
    int    time_bin_lo    ;
    int    time_bin_hi    ;

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
    TpcClusterBuilder() :
      layerGeom {nullptr}, layer {0}, side {0}, neff_electrons {0}, 
      phi_integral {0.}, time_integral {0.},
      phi_bin_lo {0}, phi_bin_hi {0}, time_bin_lo {0}, time_bin_hi {0}
    {};


    TpcClusterBuilder& operator+=(const TpcClusterBuilder& rhs);
    void reset();
    bool has_data() const;
    PairCluskeyCluster build(MapHitsetkeyUInt& cluster_cnt) const;

    ~TpcClusterBuilder(){};
};

#endif  //TRACKBASE_PADPLANEREADOUTSTRUCT_H
