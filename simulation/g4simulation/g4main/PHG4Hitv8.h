#ifndef __PHG4Hitv8_H__
#define __PHG4Hitv8_H__

#include "PHG4Hitv6.h"
#include <Geant4/G4Allocator.hh>

#ifndef __CINT__
class PHG4Hitv8;
extern G4Allocator<PHG4Hitv8> PHG4Hitv8Allocator;
#endif

/**
 * \file PHG4Hitv8.h
 * \brief Add three indices (j,k,l) to identify calorimeter tower this hit occured in
 * \author Nils Feege <nils.feege@stonybrook.edu>
 */

class PHG4Hitv8 : public PHG4Hitv6
{
public:
  PHG4Hitv8();
  PHG4Hitv8(PHG4Hit const &g4hit);

  int get_index_j() {return _idx_j;}
  int get_index_k() {return _idx_k;}
  int get_index_l() {return _idx_l;}

  void set_index_j(const int i) { _idx_j = i;}
  void set_index_k(const int i) { _idx_k = i;}
  void set_index_l(const int i) { _idx_l = i;}

  virtual void print() const;

protected:
  int _idx_j;
  int _idx_k;
  int _idx_l;

  ClassDef(PHG4Hitv8,1)
};

#endif
