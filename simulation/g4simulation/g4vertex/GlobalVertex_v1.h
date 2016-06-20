#ifndef __GLOBALVERTEX_V1_H__
#define __GLOBALVERTEX_V1_H__

#include "GlobalVertex.h"

#include <phool/PHObject.h>
#include <vector>
#include <map>
#include <iostream>

class GlobalVertex_v1 : public GlobalVertex {

public:

  GlobalVertex_v1();
  virtual ~GlobalVertex_v1();

  // PHObject virtual overloads
  
  void          identify(std::ostream& os = std::cout) const;
  void          Reset() {*this = GlobalVertex_v1();}
  int           isValid() const;
  GlobalVertex* Clone() const {return new GlobalVertex_v1(*this);}

  // vertex info
  
  unsigned int get_id() const                        {return _id;}
  void         set_id(unsigned int id)               {_id = id;}
  
  float        get_t() const                         {return _t;}
  void         set_t(float t)                        {_t = t;}

  float        get_t_err() const                     {return _t_err;}
  void         set_t_err(float t_err)                {_t_err = t_err;}
  
  float        get_x() const                         {return _pos[0];}
  void         set_x(float x)                        {_pos[0] = x;}
  
  float        get_y() const                         {return _pos[1];}
  void         set_y(float y)                        {_pos[1] = y;}

  float        get_z() const                         {return _pos[2];}
  void         set_z(float z)                        {_pos[2] = z;}

  float        get_chisq() const                     {return _chisq;}
  void         set_chisq(float chisq)                {_chisq = chisq;}

  unsigned int get_ndof() const                      {return _ndof;}
  void         set_ndof(unsigned int ndof)           {_ndof = ndof;}
  
  float        get_position(unsigned int coor) const          {return _pos[coor];}
  void         set_position(unsigned int coor, float xi)      {_pos[coor] = xi;}

  float        get_error(unsigned int i, unsigned int j) const;        //< get vertex error covar
  void         set_error(unsigned int i, unsigned int j, float value); //< set vertex error covar

  //
  // associated vertex ids methods
  //

  bool   empty_vtxids()                           const {return _vtx_ids.empty();}
  size_t size_vtxids()                            const {return _vtx_ids.size();}
  size_t count_vtxids(GlobalVertex::VTXTYPE type) const {return _vtx_ids.count(type);}

  void   clear_vtxids()                                                        {_vtx_ids.clear();}
  void   insert_vtxids(GlobalVertex::VTXTYPE type, unsigned int vtxid)         {_vtx_ids.insert(std::make_pair(type,vtxid));}
  size_t erase_vtxids(GlobalVertex::VTXTYPE type)                              {return _vtx_ids.erase(type);}
  void   erase_vtxids(GlobalVertex::VtxIter iter)                              {_vtx_ids.erase(iter);}
  void   erase_vtxids(GlobalVertex::VtxIter first, GlobalVertex::VtxIter last) {_vtx_ids.erase(first,last);}

  GlobalVertex::ConstVtxIter begin_vtxids()                          const {return _vtx_ids.begin();}
  GlobalVertex::ConstVtxIter find_vtxids(GlobalVertex::VTXTYPE type) const {return _vtx_ids.find(type);}
  GlobalVertex::ConstVtxIter end_vtxids()                            const {return _vtx_ids.end();}

  GlobalVertex::VtxIter begin_vtxids()                          {return _vtx_ids.begin();}
  GlobalVertex::VtxIter find_vtxids(GlobalVertex::VTXTYPE type) {return _vtx_ids.find(type);}
  GlobalVertex::VtxIter end_vtxids()                            {return _vtx_ids.end();}
    
private:

  unsigned int covar_index(unsigned int i, unsigned int j) const;
  
  unsigned int                                 _id;      //< unique identifier within container
  float                                        _t;       //< collision time
  float                                        _t_err;   //< collision time uncertainty
  float                                        _pos[3];  //< collision position x,y,z
  float                                        _chisq;   //< vertex fit chisq 
  unsigned int                                 _ndof;    //< degrees of freedom
  float                                        _err[6];  //< error covariance matrix (+/- cm^2)
  std::map<GlobalVertex::VTXTYPE,unsigned int> _vtx_ids; //< list of vtx ids
  
  ClassDef(GlobalVertex_v1, 1);
};

#endif

