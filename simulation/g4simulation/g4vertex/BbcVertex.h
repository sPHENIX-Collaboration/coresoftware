#ifndef __BBCVERTEX_H__
#define __BBCVERTEX_H__

#include <phool/PHObject.h>
#include <vector>
#include <set>
#include <iostream>

class BbcVertex : public PHObject {

public:
  
  BbcVertex();
  virtual ~BbcVertex();

  // PHObject virtual overloads
  
  void         identify(std::ostream& os = std::cout) const;
  void         Reset();
  int          IsValid() const;

  // vertex info
  
  unsigned int get_id() const                        {return _id;}
  void         set_id(unsigned int id)               {_id = id;}
  
  float        get_t0() const                        {return _t0;}
  void         set_t0(float t0)                      {_t0 = t0;}

  float        get_t0_err() const                    {return _t0_err;}
  void         set_t0_err(float t0_err)              {_t0_err = t0_err;}
  
  float        get_z() const                         {return _z;}
  void         set_z(float z)                        {_z = z;}

  float        get_z_err() const                     {return _z_err;}
  void         set_z_err(float z_err)                {_z_err = z_err;}
  
private:
  
  unsigned int                     _id;        //< unique identifier within container
  float                            _t0;        //< collision time
  float                            _t0_err;    //< collision time uncertainty
  float                            _z;         //< collision position z
  float                            _z_err;     //< collision position z uncertainty
  
  ClassDef(BbcVertex, 1);
};

#endif

