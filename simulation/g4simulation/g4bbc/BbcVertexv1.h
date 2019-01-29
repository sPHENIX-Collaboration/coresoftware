#ifndef G4BBC_BBCVERTEXV1_H
#define G4BBC_BBCVERTEXV1_H

#include "BbcVertex.h"

#include <phool/PHObject.h>

#include <vector>
#include <set>
#include <iostream>

class BbcVertexv1 : public BbcVertex {

public:
  
  BbcVertexv1();
  virtual ~BbcVertexv1();

  // PHObject virtual overloads
  
  void         identify(std::ostream& os = std::cout) const;
  void         Reset() {*this = BbcVertexv1();}
  int          isValid() const;
  BbcVertex*   Clone() {return new BbcVertexv1(*this);}
  
  // vertex info
  
  unsigned int get_id() const                        {return _id;}
  void         set_id(unsigned int id)               {_id = id;}
  
  float        get_t() const                         {return _t;}
  void         set_t(float t)                        {_t = t;}

  float        get_t_err() const                     {return _t_err;}
  void         set_t_err(float t_err)                {_t_err = t_err;}
  
  float        get_z() const                         {return _z;}
  void         set_z(float z)                        {_z = z;}

  float        get_z_err() const                     {return _z_err;}
  void         set_z_err(float z_err)                {_z_err = z_err;}
  
private:
  
  unsigned int                     _id;        //< unique identifier within container
  float                            _t;         //< collision time
  float                            _t_err;     //< collision time uncertainty
  float                            _z;         //< collision position z
  float                            _z_err;     //< collision position z uncertainty
  
  ClassDef(BbcVertexv1, 1);
};

#endif

