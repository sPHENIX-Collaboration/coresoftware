#ifndef __SVTXHIT_H__
#define __SVTXHIT_H__

#include <phool/PHObject.h>
#include <iostream>

class SvtxHit : public PHObject {

public:
  
  SvtxHit();
  virtual ~SvtxHit() {}

  // PHObject virtual overloads
  
  void         identify(std::ostream& os = std::cout) const;
  void         Reset() {*this = SvtxHit();}
  int          IsValid() const;

  // digitized hit info
  
  unsigned int get_id() const                        {return _id;}
  void         set_id(unsigned int id)               {_id = id;}
  
  unsigned int get_layer() const                     {return _layer;}
  void         set_layer(unsigned int layer)         {_layer = layer;}

  unsigned int get_adc() const                       {return _adc;}
  void         set_adc(unsigned int adc)             {_adc = adc;}

  float        get_e() const                         {return _e;}
  void         set_e(float e)                        {_e = e;}
  
  unsigned int get_cellid() const                    {return _cellid;}
  void         set_cellid(unsigned int cellid)       {_cellid = cellid;}
  
private:
  
  unsigned int _id;                //< unique identifier within container
  unsigned int _layer;             //< detector layer id
  unsigned int _adc;               //< digitized adc value
  float        _e;                 //< digitized energy value
  unsigned int _cellid;            //< geant4 cell object
  
  ClassDef(SvtxHit, 1);
};

#endif
