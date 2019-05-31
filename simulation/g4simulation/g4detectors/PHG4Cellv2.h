#ifndef G4DETECTORS_PHG4CELLV2_H
#define G4DETECTORS_PHG4CELLV2_H

#include "PHG4Cell.h"
#include "PHG4CellDefs.h"

#include <g4main/PHG4HitDefs.h>  // for keytype

#include <iostream>
#include <map>
#include <utility>               // for make_pair

//! specialized cells for TPC operations
class PHG4Cellv2 : public PHG4Cell
{
 public:
  PHG4Cellv2();
  explicit PHG4Cellv2(const PHG4CellDefs::keytype g4cellid);
  virtual ~PHG4Cellv2();

  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset();

  void set_cellid(const PHG4CellDefs::keytype i) { cellid = i; }
  PHG4CellDefs::keytype get_cellid() const { return cellid; }
  bool has_binning(const PHG4CellDefs::CellBinning binning) const;
  short int get_detid() const;

  void add_edep(const PHG4HitDefs::keytype g4hitid, const float edep);

  EdepConstRange get_g4hits()
  {
    return std::make_pair(hitedeps.begin(), hitedeps.end());
  }

  void add_edep(const float f) { _edep += f; }
  double get_edep() const { return _edep; }
  //  tpctod* get_train_of_digits() {return &trainOfDigits;}

  void print() const;

 protected:
  PHG4CellDefs::keytype cellid;
  EdepMap hitedeps;

  float _edep;
  //  tpctod trainOfDigits;

  ClassDef(PHG4Cellv2, 1)
};

#endif
