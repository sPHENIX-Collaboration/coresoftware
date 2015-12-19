#ifndef __PHG4SHOWER_V1_H__
#define __PHG4SHOWER_V1_H__

#include "PHG4Shower.h"

#include "PHG4HitDefs.h"

#include <phool/PHObject.h>
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
  
  unsigned int get_id() const                 {return _id;}
  void         set_id(unsigned int id)        {_id = id;}

  int          get_primary_id() const         {return _primary_id;}
  void         set_primary_id(int primary_id) {_primary_id = primary_id;}

  int          get_parent_shower_id() const   {return _parent_shower_id;}
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
  
  float        get_edep(int volume) const;
  void         set_edep(int volume, float edep) {_edep[volume] = edep;}

  float        get_eion(int volume) const;
  void         set_eion(int volume, float eion) {_eion[volume] = eion;}

  float        get_light_yield(int volume) const;
  void         set_light_yield(int volume, float light_yield) {_light_yield[volume] = light_yield;}

  // container methods for ids
  void                            add_g4particle_id(int id)  {_g4particle_ids.insert(id);}
  PHG4Shower::ParticleIdIter      begin_g4particle_id() {return _g4particle_ids.begin();}
  PHG4Shower::ParticleIdConstIter begin_g4particle_id() const {return _g4particle_ids.begin();}
  PHG4Shower::ParticleIdIter      end_g4particle_id() {return _g4particle_ids.end();}
  PHG4Shower::ParticleIdConstIter end_g4particle_id() const {return _g4particle_ids.end();}
  size_t                          remove_g4particle_id(int id) {return _g4particle_ids.erase(id);}
  
  void                       add_g4hit_id(int volume,PHG4HitDefs::keytype id) {_g4hit_ids[volume].insert(id);}
  PHG4Shower::HitIdIter      begin_g4hit_id() {return _g4hit_ids.begin();}
  PHG4Shower::HitIdConstIter begin_g4hit_id() const {return _g4hit_ids.begin();}
  PHG4Shower::HitIdIter      end_g4hit_id() {return _g4hit_ids.end();}
  PHG4Shower::HitIdConstIter end_g4hit_id() const {return _g4hit_ids.end();}
  size_t                     remove_g4hit_id(int volume,int id) {return _g4hit_ids[volume].erase(id);}

private:
  
  unsigned int covar_index(unsigned int i, unsigned int j) const;
  
  unsigned int         _id;               //< unique identifier within container
  int                  _primary_id;       //< association of shower to primary particle id
  int                  _parent_shower_id; //< association of shower to parent shower id if present
  float                _pos[3];           //< mean position of the shower
  float                _covar[6];         //< covariance of shower positions
  std::map<int, float> _edep;             //< energy deposit in different volumes
  std::map<int, float> _eion;             //< ionization energy in different volumes
  std::map<int, float> _light_yield;      //< light yield in different volumes

  std::set<int> _g4particle_ids;
  std::map<int,std::set<PHG4HitDefs::keytype> > _g4hit_ids;
  
  ClassDef(PHG4Showerv1, 1);
};

#endif

