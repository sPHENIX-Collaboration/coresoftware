/*!
 * \file Jetv1.h
 * \brief Versionize the Jet object that make by Mike McCumber
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4JET_JETV1_H
#define G4JET_JETV1_H

#include "Jet.h"

/*!
 * \brief Jetv1
 */
class Jetv1 : public Jet {
  
public:
  
  Jetv1();
  virtual ~Jetv1() {}

  // PHObject virtual overloads

  void         identify(std::ostream& os = std::cout) const;
  void         Reset();
  int          isValid() const;
  Jet*         Clone() const;

  // jet info

  unsigned int get_id() const            {return _id;}
  void         set_id(unsigned int id)   {_id = id;}

  float        get_px() const            {return _mom[0];}
  void         set_px(float px)          {_mom[0] = px;}

  float        get_py() const            {return _mom[1];}
  void         set_py(float py)          {_mom[1] = py;}

  float        get_pz() const            {return _mom[2];}
  void         set_pz(float pz)          {_mom[2] = pz;}

  float        get_e() const             {return _e;}
  void         set_e(float e)            {_e = e;}

  float        get_p() const;
  float        get_pt() const;
  float        get_et() const;
  float        get_eta() const;
  float        get_phi() const;
  float        get_mass() const;
  float        get_mass2() const;
  
  // extended jet info
  
  bool  has_property(Jet::PROPERTY prop_id) const;
  float get_property(Jet::PROPERTY prop_id) const;
  void  set_property(Jet::PROPERTY prop_id, float value);
  void  print_property(std::ostream& os) const;
  
  //
  // clustered component methods (multimap interface based)
  // source type id --> unique id within that storage
  //
  bool      empty_comp() const           {return _comp_ids.empty();}
  size_t    size_comp() const            {return _comp_ids.size();}
  size_t    count_comp(SRC source) const {return _comp_ids.count(source);}

  void      clear_comp()                                {_comp_ids.clear();}
  void      insert_comp(SRC source,unsigned int compid) {_comp_ids.insert(std::make_pair(source,compid));}
  size_t    erase_comp(SRC source)                      {return _comp_ids.erase(source);}
  void      erase_comp(Iter iter)                       {_comp_ids.erase(iter); return;}
  void      erase_comp(Iter first, Iter last)           {_comp_ids.erase(first,last); return;}

  ConstIter begin_comp() const                 {return _comp_ids.begin();}
  ConstIter lower_bound_comp(SRC source) const {return _comp_ids.lower_bound(source);}
  ConstIter upper_bound_comp(SRC source) const {return _comp_ids.upper_bound(source);}
  ConstIter find(SRC source) const             {return _comp_ids.find(source);}
  ConstIter end_comp() const                   {return _comp_ids.end();}

  Iter begin_comp()                 {return _comp_ids.begin();}
  Iter lower_bound_comp(SRC source) {return _comp_ids.lower_bound(source);}
  Iter upper_bound_comp(SRC source) {return _comp_ids.upper_bound(source);}
  Iter find(SRC source)             {return _comp_ids.find(source);}
  Iter end_comp()                   {return _comp_ids.end();}

private:

  /// unique identifier within container
  unsigned int _id;

  /// jet momentum vector (px,py,pz)
  float _mom[3];

  /// jet energy
  float _e;

  /// source id -> component id
  typ_comp_ids _comp_ids;

  typedef std::map<Jet::PROPERTY, float> typ_property_map;
  /// map that contains extra properties
  typ_property_map _property_map;

  ClassDef(Jetv1, 1);
};

#endif // G4JET_JETV1_H
