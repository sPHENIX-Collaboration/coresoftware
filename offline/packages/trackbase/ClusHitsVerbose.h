#ifndef CLUSHITSVERBOSE__H
#define CLUSHITSVERBOSE__H

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <map>
#include <vector>
#include <array>

class ClusHitsVerbose : public PHObject
{
 public:
  using BinData = std::pair<int,int>;  // index + energy for a given bin (MVTX, INTT, or TPC)
  using Vector  = std::vector<BinData>;// listing in phi and z (will need two of them)
  using Map     = std::map<TrkrDefs::cluskey, std::array<Vector,4>>;
  // ^ what is stored:: first: Phi data,  Z data, cut Phi data, cut Z data */
  void Reset() override {}

  virtual bool hasClusKey (TrkrDefs::cluskey) const { return false; };
  virtual Vector& phiBins    (TrkrDefs::cluskey);
  virtual Vector& zBins      (TrkrDefs::cluskey);
  virtual Vector& phiCutBins (TrkrDefs::cluskey);
  virtual Vector& zCutBins   (TrkrDefs::cluskey);
  virtual Map& getMap();

  // convenience classes; 
  // ROOT's default library has vector<vector<int>> but not vec<vec<pair<int,int>>>
  using VecInt   = std::vector<int>;
  using PairVector = std::pair<VecInt,VecInt>;//
  virtual PairVector phiBins_pvecIE(TrkrDefs::cluskey); // IE for integer-energy
  virtual PairVector phiCutBins_pvecIE(TrkrDefs::cluskey);
  virtual PairVector zBins_pvecIE(TrkrDefs::cluskey);
  virtual PairVector zCutBins_pvecIE(TrkrDefs::cluskey);
  virtual void addPhiHit(int, int) {return;}
  virtual void addZHit(int, int) {return;}
  virtual void addPhiCutHit(int, int) {return;}
  virtual void addZCutHit(int, int) {return;}
  virtual void push_hits (TrkrDefs::cluskey) {return;}

  // PHObject virtual overload
  void identify(std::ostream& os = std::cout) const override
  {
    os << "ClusHitsVerbose base class" << std::endl;
  };

 protected:
  ClusHitsVerbose() = default;
  ClassDefOverride(ClusHitsVerbose, 1)
};

#endif  // G4TRACKING_CLUSHITSVERBOSE_H
