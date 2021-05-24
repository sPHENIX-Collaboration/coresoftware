// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4SHOWERV1_H
#define G4MAIN_PHG4SHOWERV1_H

#include "PHG4Shower.h"

#include "PHG4HitDefs.h"

#include <cstddef>       // for size_t
#include <iostream>
#include <map>
#include <set>

class PHG4Showerv1 : public PHG4Shower
{
 public:
  PHG4Showerv1();
  ~PHG4Showerv1() override {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override;
  PHG4Shower* CloneMe() const override { return (new PHG4Showerv1(*this)); }
  void Reset() override { *this = PHG4Showerv1(); }
  int isValid() const override;

  // shower info

  int get_id() const override { return _id; }
  void set_id(int id) override { _id = id; }

  int get_parent_particle_id() const override { return _parent_particle_id; }
  void set_parent_particle_id(int parent_particle_id) override { _parent_particle_id = parent_particle_id; }

  int get_parent_shower_id() const override { return _parent_shower_id; }
  void set_parent_shower_id(int parent_shower_id) override { _parent_shower_id = parent_shower_id; }

  float get_x() const override { return _pos[0]; }
  void set_x(float x) override { _pos[0] = x; }

  float get_y() const override { return _pos[1]; }
  void set_y(float y) override { _pos[1] = y; }

  float get_z() const override { return _pos[2]; }
  void set_z(float z) override { _pos[2] = z; }

  float get_position(unsigned int coor) const override { return _pos[coor]; }
  void set_position(unsigned int coor, float xi) override { _pos[coor] = xi; }

  float get_covar(unsigned int i, unsigned int j) const override;
  void set_covar(unsigned int i, unsigned int j, float entry) override;

  unsigned int get_nhits(int volume) const override;
  void set_nhits(int volume, unsigned int nhits) override { _nhits[volume] = nhits; }

  double get_edep() const override;
  float get_edep(int volume) const override;
  void set_edep(int volume, float edep) override { _edep[volume] = edep; }

  double get_eion() const override;
  float get_eion(int volume) const override;
  void set_eion(int volume, float eion) override { _eion[volume] = eion; }

  float get_light_yield(int volume) const override;
  void set_light_yield(int volume, float light_yield) override { _light_yield[volume] = light_yield; }

  float get_eh_ratio(int volume) const override;
  void set_eh_ratio(int volume, float eh_ratio) override { _eh_ratio[volume] = eh_ratio; }

  // container methods for ids
  bool empty_g4particle_id() const override { return _g4particle_ids.empty(); }
  size_t size_g4particle_id() const override { return _g4particle_ids.size(); }
  void add_g4particle_id(int id) override { _g4particle_ids.insert(id); }
  PHG4Shower::ParticleIdIter begin_g4particle_id() override { return _g4particle_ids.begin(); }
  PHG4Shower::ParticleIdConstIter begin_g4particle_id() const override { return _g4particle_ids.begin(); }
  PHG4Shower::ParticleIdIter end_g4particle_id() override { return _g4particle_ids.end(); }
  PHG4Shower::ParticleIdConstIter end_g4particle_id() const override { return _g4particle_ids.end(); }
  size_t remove_g4particle_id(int id) override { return _g4particle_ids.erase(id); }
  void clear_g4particle_id() override { return _g4particle_ids.clear(); }
  const ParticleIdSet& g4particle_ids() const override { return _g4particle_ids; }

  bool empty_g4vertex_id() const override { return _g4vertex_ids.empty(); }
  size_t size_g4vertex_id() const override { return _g4vertex_ids.size(); }
  void add_g4vertex_id(int id) override { _g4vertex_ids.insert(id); }
  PHG4Shower::VertexIdIter begin_g4vertex_id() override { return _g4vertex_ids.begin(); }
  PHG4Shower::VertexIdConstIter begin_g4vertex_id() const override { return _g4vertex_ids.begin(); }
  PHG4Shower::VertexIdIter end_g4vertex_id() override { return _g4vertex_ids.end(); }
  PHG4Shower::VertexIdConstIter end_g4vertex_id() const override { return _g4vertex_ids.end(); }
  size_t remove_g4vertex_id(int id) override { return _g4vertex_ids.erase(id); }
  void clear_g4vertex_id() override { return _g4vertex_ids.clear(); }
  const VertexIdSet& g4vertex_ids() const override { return _g4vertex_ids; }

  bool empty_g4hit_id() const override { return _g4hit_ids.empty(); }
  size_t size_g4hit_id() const override { return _g4hit_ids.size(); }
  void add_g4hit_id(int volume, PHG4HitDefs::keytype id) override { _g4hit_ids[volume].insert(id); }
  PHG4Shower::HitIdIter begin_g4hit_id() override { return _g4hit_ids.begin(); }
  PHG4Shower::HitIdConstIter begin_g4hit_id() const override { return _g4hit_ids.begin(); }
  PHG4Shower::HitIdIter find_g4hit_id(int volume) override { return _g4hit_ids.find(volume); }
  PHG4Shower::HitIdConstIter find_g4hit_id(int volume) const override { return _g4hit_ids.find(volume); }
  PHG4Shower::HitIdIter end_g4hit_id() override { return _g4hit_ids.end(); }
  PHG4Shower::HitIdConstIter end_g4hit_id() const override { return _g4hit_ids.end(); }
  size_t remove_g4hit_id(int volume, PHG4HitDefs::keytype id) override { return _g4hit_ids[volume].erase(id); }
  size_t remove_g4hit_volume(int volume) override { return _g4hit_ids.erase(volume); }
  void clear_g4hit_id() override { return _g4hit_ids.clear(); }
  const HitIdMap& g4hit_ids() const override { return _g4hit_ids; }

 private:
  unsigned int covar_index(unsigned int i, unsigned int j) const;

  int _id;                             //< unique identifier within container
  int _parent_particle_id;             //< association of shower to parent particle id
  int _parent_shower_id;               //< association of shower to parent shower id
  float _pos[3];                       //< mean position of the shower hits
  float _covar[6];                     //< covariance of shower hits
  std::map<int, unsigned int> _nhits;  //< number of hits in different volumes
  std::map<int, float> _edep;          //< energy deposit in different volumes
  std::map<int, float> _eion;          //< ionization energy in different volumes
  std::map<int, float> _light_yield;   //< light yield in different volumes
  std::map<int, float> _eh_ratio;      //< electron/hadron ratio of energy in different volumes

  // these containers are cleared during dst reduction, but are available in full dsts
  ParticleIdSet _g4particle_ids;       //< contained secondary particle ids
  VertexIdSet _g4vertex_ids;           //< contained secondary vertex ids
  HitIdMap _g4hit_ids;                 //< contained hit ids

  ClassDefOverride(PHG4Showerv1, 1);
};

#endif
