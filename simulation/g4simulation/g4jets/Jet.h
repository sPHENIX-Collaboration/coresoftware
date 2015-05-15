#ifndef __JET_H__
#define __JET_H__

#include <phool/PHObject.h>
#include <map>
#include <iostream>

class Jet : public PHObject {

public:

  enum ALGO {NONE,ANTIKT,KT,CONE};
  
  enum SRC {TRACKS,
	    EM_TOWERS,EM_CLUSTERS,
	    HCALIN_TOWERS,HCALIN_CLUSTERS,
	    HCALOUT_TOWERS,HCALOUT_CLUSTERS};
  
  typedef std::multimap< SRC, unsigned int >::const_iterator ConstIter;
  typedef std::multimap< SRC, unsigned int >::iterator       Iter; 
  
  Jet(); 
  virtual ~Jet() {}

  // PHObject virtual overloads
  
  void         identify(std::ostream& os = std::cout) const;
  void         Reset();
  int          IsValid() const;

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
   
private:
  
  unsigned int _id;                //< unique identifier within container
  float _mom[3];                   //< jet momentum vector (px,py,pz)
  float _e;                        //< jet energy

  std::multimap< SRC, unsigned int > _comp_ids; //< source id -> component id 
  
  ClassDef(Jet, 1);
};

#endif
