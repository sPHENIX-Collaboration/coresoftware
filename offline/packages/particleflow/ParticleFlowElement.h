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

class SvtxTrack;
class RawCluster;

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
  ~ParticleFlowElement() override {}

  void identify(std::ostream& os = std::cout) const override;
  int isValid() const override { return 0; }

  virtual unsigned int get_id() const { return 0xFFFFFFFF; }
  virtual void set_id(unsigned int) { return; }

  virtual ParticleFlowElement::PFLOWTYPE get_type() const {return ParticleFlowElement::PFLOWTYPE::UNASSIGNED; }
  virtual void set_type( ParticleFlowElement::PFLOWTYPE ) { return; }

  virtual float get_px() const { return NAN; }
  virtual void set_px(float) { return; }

  virtual float get_py() const { return NAN; }
  virtual void set_py(float) { return; }

  virtual float get_pz() const { return NAN; }
  virtual void set_pz(float) { return; }

  virtual float get_e() const { return NAN; }
  virtual void set_e(float) { return; }

  virtual SvtxTrack* get_track() const { return nullptr; }
  virtual void set_track(SvtxTrack*) { return; }
  
  virtual std::vector<RawCluster*> get_eclusters() const { return std::vector<RawCluster*>(); }
  virtual void set_eclusters(const std::vector<RawCluster*>&) { return; }
  
  virtual RawCluster* get_hcluster() const { return nullptr; }
  virtual void set_hcluster(RawCluster*) { return; }

  virtual float get_p() const { return NAN; }
  virtual float get_pt() const { return NAN; }
  virtual float get_et() const { return NAN; }
  virtual float get_eta() const { return NAN; }
  virtual float get_phi() const { return NAN; }
  virtual float get_mass() const { return NAN; }

  ClassDefOverride(ParticleFlowElement, 1);
};

#endif
