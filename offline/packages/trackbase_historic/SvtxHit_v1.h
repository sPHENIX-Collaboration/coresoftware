#ifndef TRACKBASEHISTORIC_SVTXHITV1_H
#define TRACKBASEHISTORIC_SVTXHITV1_H

#include "SvtxHit.h"

#include <g4detectors/PHG4CellDefs.h>

#include <iostream>

class SvtxHit_v1 : public SvtxHit
{
 public:
  SvtxHit_v1();
  virtual ~SvtxHit_v1() {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const;
  void Reset() { *this = SvtxHit_v1(); }
  int isValid() const;
  SvtxHit* Clone() const { return (new SvtxHit_v1(*this)); }

  // digitized hit info

  unsigned int get_id() const { return _id; }
  void set_id(unsigned int id) { _id = id; }

  unsigned int get_layer() const { return _layer; }
  void set_layer(unsigned int layer) { _layer = layer; }

  unsigned int get_adc() const { return _adc; }
  void set_adc(unsigned int adc) { _adc = adc; }

  float get_e() const { return _e; }
  void set_e(float e) { _e = e; }

  PHG4CellDefs::keytype get_cellid() const { return _cellid; }
  void set_cellid(PHG4CellDefs::keytype cellid) { _cellid = cellid; }

 private:
  unsigned int _id;               //< unique identifier within container
  unsigned int _layer;            //< detector layer id
  unsigned int _adc;              //< digitized adc value
  float _e;                       //< digitized energy value
  PHG4CellDefs::keytype _cellid;  //< geant4 cell object

  ClassDef(SvtxHit_v1, 2);
};

#endif
