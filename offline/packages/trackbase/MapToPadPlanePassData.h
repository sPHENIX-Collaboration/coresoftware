#ifndef TRACKBASE_MAPTOPADPLANEPASSDATA_H
#define TRACKBASE_MAPTOPADPLANEPASSDATA_H

#include <phool/PHObject.h>
#include "TrkrDefs.h"
  // This is really just a structure used in simulation/g4simulation/g4tpc/PHG4TpcElectronDrift.cc
  // to collect necessary statistics from simulation/g4simulation/g4tpc/PHG4TpcPadPlaneReadout.cc
  // 
  // It could have been just std::pair< std::array<int,6>, std::pair<double,double>>, but 
  // having the structure to name those subtypes.


class MapToPadPlanePassData: public PHObject
{
  public:
    TrkrDefs::hitsetkey hitsetkey;
    int    neff_electrons ;
    double phi_integral   ;
    double time_integral  ;
    int    phi_bin_lo     ;
    int    phi_bin_hi     ;
    int    time_bin_lo    ;
    int    time_bin_hi    ;

    MapToPadPlanePassData( TrkrDefs::hitsetkey _hitsetkey, int _neff_electrons, double _phi_integral, double _time_integral,
        int _phi_bin_lo,  int _phi_bin_hi,
        int _time_bin_lo, int _time_bin_hi) :
          hitsetkey          { _hitsetkey          },
          neff_electrons { _neff_electrons },
          phi_integral   { _phi_integral   },
          time_integral  { _time_integral  },
          phi_bin_lo     { _phi_bin_lo     },
          phi_bin_hi     { _phi_bin_hi     },
          time_bin_lo    { _time_bin_lo    },
          time_bin_hi    { _time_bin_hi    }
    {};
    MapToPadPlanePassData() :
      hitsetkey {0}, neff_electrons {0}, phi_integral {0.}, time_integral {0.},
      phi_bin_lo {0}, phi_bin_hi {0}, time_bin_lo {0}, time_bin_hi {0}
    {};


    MapToPadPlanePassData& operator+=(const MapToPadPlanePassData& rhs);
    void reset();

    ~MapToPadPlanePassData() = default;
  /* protected: */

  private:
    ClassDefOverride(MapToPadPlanePassData, 1);
};



#endif  //TRACKBASE_PADPLANEREADOUTSTRUCT_H
