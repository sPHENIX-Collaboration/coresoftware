#ifndef G4MAIN_PHG4SHOWERV1_H
#define G4MAIN_PHG4SHOWERV1_H

#include "PHG4Shower.h"

#include "PHG4HitDefs.h"

#include <cstddef>
#include <map>
#include <set>
#include <iostream>

class PHG4Showerv1 : public PHG4Shower {

public:
  
  PHG4Showerv1();
  virtual ~PHG4Showerv1() {}

  // PHObject virtual overloads
   
  void         identify(std::ostream& os = std::cout) const;
  PHG4Shower*  Clone() const {return (new PHG4Showerv1(*this));}
  void         Reset() {*this = PHG4Showerv1();}
  int          isValid() const;

  // shower info
  
  int          get_id() const {return _id;}
  void         set_id(int id) {_id = id;}

  int          get_parent_particle_id() const {return _parent_particle_id;}
  void         set_parent_particle_id(int parent_particle_id) {_parent_particle_id = parent_particle_id;}
  
  int          get_parent_shower_id() const {return _parent_shower_id;}
  void         set_parent_shower_id(int parent_shower_id) {_parent_shower_id = parent_shower_id;}
  
  float        get_x() const                  {return _pos[0];}
  void         set_x(float x)                 {_pos[0] = x;}

  float        get_y() const                  {return _pos[1];}
  void         set_y(float y)                 {_pos[1] = y;}

  float        get_z() const                  {return _pos[2];}
  void         set_z(float z)                 {_pos[2] = z;}

  float        get_position(unsigned int coor) const     {return _pos[coor];}
  void         set_position(unsigned int coor, float xi) {_pos[coor] = xi;}
  
  float        get_covar(unsigned int i, unsigned int j) const;
  void         set_covar(unsigned int i, unsigned int j, float entry);

  unsigned int get_nhits(int volume) const;
  void         set_nhits(int volume, unsigned int nhits) {_nhits[volume] = nhits;}
  
  double       get_edep() const ;
  float        get_edep(int volume) const;
  void         set_edep(int volume, float edep) {_edep[volume] = edep;}

  double       get_eion() const ;
  float        get_eion(int volume) const;
  void         set_eion(int volume, float eion) {_eion[volume] = eion;}

  float        get_light_yield(int volume) const;
  void         set_light_yield(int volume, float light_yield) {_light_yield[volume] = light_yield;}

  float        get_eh_ratio(int volume) const;
  void         set_eh_ratio(int volume, float eh_ratio) {_eh_ratio[volume] = eh_ratio;}
  
  // container methods for ids
  bool                            empty_g4particle_id() const {return _g4particle_ids.empty();}
  size_t                          size_g4particle_id() const {return _g4particle_ids.size();}
  void                            add_g4particle_id(int id)  {_g4particle_ids.insert(id);}
  PHG4Shower::ParticleIdIter      begin_g4particle_id() {return _g4particle_ids.begin();}
  PHG4Shower::ParticleIdConstIter begin_g4particle_id() const {return _g4particle_ids.begin();}
  PHG4Shower::ParticleIdIter      end_g4particle_id() {return _g4particle_ids.end();}
  PHG4Shower::ParticleIdConstIter end_g4particle_id() const {return _g4particle_ids.end();}
  size_t                          remove_g4particle_id(int id) {return _g4particle_ids.erase(id);}
  void                            clear_g4particle_id() {return _g4particle_ids.clear();}

  bool                          empty_g4vertex_id() const {return _g4vertex_ids.empty();}
  size_t                        size_g4vertex_id() const {return _g4vertex_ids.size();}
  void                          add_g4vertex_id(int id)  {_g4vertex_ids.insert(id);}
  PHG4Shower::VertexIdIter      begin_g4vertex_id() {return _g4vertex_ids.begin();}
  PHG4Shower::VertexIdConstIter begin_g4vertex_id() const {return _g4vertex_ids.begin();}
  PHG4Shower::VertexIdIter      end_g4vertex_id() {return _g4vertex_ids.end();}
  PHG4Shower::VertexIdConstIter end_g4vertex_id() const {return _g4vertex_ids.end();}
  size_t                        remove_g4vertex_id(int id) {return _g4vertex_ids.erase(id);}
  void                          clear_g4vertex_id() {return _g4vertex_ids.clear();}
  
  bool                       empty_g4hit_id() const {return _g4hit_ids.empty();}
  size_t                     size_g4hit_id() const {return _g4hit_ids.size();}
  void                       add_g4hit_id(int volume,PHG4HitDefs::keytype id) {_g4hit_ids[volume].insert(id);}
  PHG4Shower::HitIdIter      begin_g4hit_id() {return _g4hit_ids.begin();}
  PHG4Shower::HitIdConstIter begin_g4hit_id() const {return _g4hit_ids.begin();}
  PHG4Shower::HitIdIter      find_g4hit_id(int volume) {return _g4hit_ids.find(volume);}
  PHG4Shower::HitIdConstIter find_g4hit_id(int volume) const {return _g4hit_ids.find(volume);}
  PHG4Shower::HitIdIter      end_g4hit_id() {return _g4hit_ids.end();}
  PHG4Shower::HitIdConstIter end_g4hit_id() const {return _g4hit_ids.end();}
  size_t                     remove_g4hit_id(int volume, PHG4HitDefs::keytype id) {return _g4hit_ids[volume].erase(id);}
  size_t                     remove_g4hit_volume(int volume) {return _g4hit_ids.erase(volume);}
  void                       clear_g4hit_id() {return _g4hit_ids.clear();}
  
private:
  
  unsigned int covar_index(unsigned int i, unsigned int j) const;
  
  int                  _id;                 //< unique identifier within container
  int                  _parent_particle_id; //< association of shower to parent particle id
  int                  _parent_shower_id;   //< association of shower to parent shower id
  float                _pos[3];             //< mean position of the shower hits
  float                _covar[6];           //< covariance of shower hits
  std::map<int, unsigned int> _nhits;       //< number of hits in different volumes
  std::map<int, float> _edep;               //< energy deposit in different volumes
  std::map<int, float> _eion;               //< ionization energy in different volumes
  std::map<int, float> _light_yield;        //< light yield in different volumes
  std::map<int, float> _eh_ratio;           //< electron/hadron ratio of energy in different volumes

  // these containers are cleared during dst reduction, but are available in full dsts
  std::set<int> _g4particle_ids; //< contained secondary particle ids
  std::set<int> _g4vertex_ids;   //< contained secondary vertex ids
  std::map<int,std::set<unsigned long long> > _g4hit_ids; //< contained hit ids
  
  ClassDef(PHG4Showerv1, 1);
};

#endif

