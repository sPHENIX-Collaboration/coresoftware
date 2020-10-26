#ifndef PARTICLEFLOW_PARTICLEFLOWELEMENTV1_H
#define PARTICLEFLOW_PARTICLEFLOWELEMENTV1_H

//===========================================================
/// \file ParticleFlowElementv1.h
/// \brief v1 Implementation of ParticleFlowElement base class 
/// \author Dennis V. Perepelitsa
//===========================================================

#include "TObject.h"

#include <iostream>
#include <utility>   // for pair, make_pair

#include "ParticleFlowElement.h"

class ParticleFlowElementv1 : public ParticleFlowElement
{
 public:
  ParticleFlowElementv1();
  virtual ~ParticleFlowElementv1() {}
  
  // PHObject virtual overloads
  
  void identify(std::ostream& os = std::cout) const;
  void Reset();
  int isValid() const;
  
  // pflow element info

  ParticleFlowElement::PFLOWTYPE get_type() const {return _type; }
  void set_type( ParticleFlowElement::PFLOWTYPE type ) { _type = type; }
  
  unsigned int get_id() const { return _id; }
  void set_id(unsigned int id) { _id = id; }
  
  float get_px() const { return _mom[0]; }
  void set_px(float px) { _mom[0] = px; }
  
  float get_py() const { return _mom[1]; }
  void set_py(float py) { _mom[1] = py; }
  
  float get_pz() const { return _mom[2]; }
  void set_pz(float pz) { _mom[2] = pz; }
  
  float get_e() const { return _e; }
  void set_e(float e) { _e = e; }
  
  float get_p() const;
  float get_pt() const;
  float get_et() const;
  float get_eta() const;
  float get_phi() const;
  float get_mass() const;
  
 private:
  /// unique identifier within container
  unsigned int _id;

  // particle flow type 
  ParticleFlowElement::PFLOWTYPE _type;
  
  /// pflow momentum vector (px,py,pz)
  float _mom[3];
  
  /// pflow energy
  float _e;
  
  ClassDef(ParticleFlowElementv1, 1);
};

#endif
