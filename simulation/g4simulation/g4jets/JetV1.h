// $Id: $                                                                                             

/*!
 * \file JetV1.h
 * \brief Versionize the Jet object that make by Mike McCumber
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef JETV1_H_
#define JETV1_H_

#include "Jet.h"

/*!
 * \brief JetV1
 */
class JetV1 : public Jet
{
public:
  JetV1();
  virtual
  ~JetV1();

  // PHObject virtual overloads

  void         identify(std::ostream& os = std::cout) const;
  void         Reset();
  int          isValid() const;

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
  void      erase_comp(Iter iter)                       {return _comp_ids.erase(iter);}
  void      erase_comp(Iter first, Iter last)           {return _comp_ids.erase(first,last);}

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

  //! whether a property exists
  bool  has_property(PROPERTY prop_id) const;
  //! get property value
  float   get_property(PROPERTY prop_id) const;
  //! set property value
  void    set_property(PROPERTY prop_id, float value);
  void print_property(ostream& os) const;

private:

  //! unique identifier within container
  unsigned int _id;

  //! jet momentum vector (px,py,pz)
  float _mom[3];

  //! jet energy
  float _e;

  //! source id -> component id
  typ_comp_ids _comp_ids;

  typedef std::map<PROPERTY, float> typ_property_map;

  //! map that contains extra properties
  typ_property_map _property_map;

  ClassDef(JetV1, 1);

};

#endif /* JETV1_H_ */
