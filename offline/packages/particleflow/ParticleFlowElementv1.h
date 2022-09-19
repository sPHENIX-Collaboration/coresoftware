#ifndef PARTICLEFLOW_PARTICLEFLOWELEMENTV1_H
#define PARTICLEFLOW_PARTICLEFLOWELEMENTV1_H

//===========================================================
/// \file ParticleFlowElementv1.h
/// \brief v1 Implementation of ParticleFlowElement base class 
/// \author Dennis V. Perepelitsa
//===========================================================

#include "ParticleFlowElement.h"

#include <iostream>
#include <utility>   // for pair, make_pair

class SvtxTrack;
class RawCluster;

class ParticleFlowElementv1 : public ParticleFlowElement
{
 public:
  ParticleFlowElementv1();
  ~ParticleFlowElementv1() override {}
  
  // PHObject virtual overloads
  
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  
  // pflow element info

  ParticleFlowElement::PFLOWTYPE get_type() const override {return _type; }
  void set_type( ParticleFlowElement::PFLOWTYPE type ) override { _type = type; }
  
  unsigned int get_id() const override { return _id; }
  void set_id(unsigned int id) override { _id = id; }
  
  float get_px() const override { return _mom[0]; }
  void set_px(float px) override { _mom[0] = px; }
  
  float get_py() const override { return _mom[1]; }
  void set_py(float py) override { _mom[1] = py; }
  
  float get_pz() const override { return _mom[2]; }
  void set_pz(float pz) override { _mom[2] = pz; }
  
  float get_e() const override { return _e; }
  void set_e(float e) override { _e = e; }

  SvtxTrack* get_track() const override { return _track; }
  void set_track(SvtxTrack* track) override { _track = track; }
  
  std::vector<RawCluster*> get_eclusters() const override { return _ecluster; }
  void set_eclusters(const std::vector<RawCluster*>& ecluster) override { _ecluster = ecluster; }
  
  RawCluster* get_hcluster() const override { return _hcluster; }
  void set_hcluster(RawCluster* hcluster) override { _hcluster = hcluster; }

  float get_p() const override;
  float get_pt() const override;
  float get_et() const override;
  float get_eta() const override;
  float get_phi() const override;
  float get_mass() const override;
  
 private:
  /// unique identifier within container
  unsigned int _id;

  // particle flow type 
  ParticleFlowElement::PFLOWTYPE _type;
  
  /// pflow momentum vector (px,py,pz)
  float _mom[3];
  
  /// pflow energy
  float _e;
  
  SvtxTrack* _track = nullptr;
  std::vector<RawCluster*> _ecluster;
  RawCluster* _hcluster = nullptr;
  ClassDefOverride(ParticleFlowElementv1, 1);
};

#endif
