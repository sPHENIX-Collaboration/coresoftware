#ifndef PARTICLEFLOW_PARTICLEFLOWELEMENT_H
#define PARTICLEFLOW_PARTICLEFLOWELEMENT_H

//===========================================================
/// \file ParticleFlowElement.h
/// \brief Base class for simple particle flow element objects
/// \author Dennis V. Perepelitsa
//===========================================================

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>

class ParticleFlowElement : public PHObject
{
 public:
  // enums can be extended with new values, but values not altered

  enum PFLOWTYPE
  {
    UNASSIGNED = -1,
    MATCHED_CHARGED_HADRON = 0,
    UNMATCHED_CHARGED_HADRON = 1,
    UNMATCHED_EM_PARTICLE = 2,
    UNMATCHED_NEUTRAL_HADRON = 3,
    LEFTOVER_EM_PARTICLE = 4
  };

  ParticleFlowElement() {}
  virtual ~ParticleFlowElement() {}

  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset() { return; }
  virtual int isValid() const { return 0; }

  virtual unsigned int get_id() const { return 0xFFFFFFFF; }
  virtual void set_id(unsigned int id) { return; }

  virtual ParticleFlowElement::PFLOWTYPE get_type() const {return ParticleFlowElement::PFLOWTYPE::UNASSIGNED; }
  virtual void set_type( ParticleFlowElement::PFLOWTYPE ) { return; }

  virtual float get_px() const { return NAN; }
  virtual void set_px(float px) { return; }

  virtual float get_py() const { return NAN; }
  virtual void set_py(float py) { return; }

  virtual float get_pz() const { return NAN; }
  virtual void set_pz(float pz) { return; }

  virtual float get_e() const { return NAN; }
  virtual void set_e(float e) { return; }

  virtual float get_p() const { return NAN; }
  virtual float get_pt() const { return NAN; }
  virtual float get_et() const { return NAN; }
  virtual float get_eta() const { return NAN; }
  virtual float get_phi() const { return NAN; }
  virtual float get_mass() const { return NAN; }

  ClassDef(ParticleFlowElement, 1);
};

#endif
