#ifndef __SVTXHIT_H__
#define __SVTXHIT_H__

#include <TObject.h>
#include <g4detectors/PHG4CellDefs.h>
#include <iostream>
#include <limits.h>

class SvtxHit : public TObject {

public:
  
  virtual ~SvtxHit() {}

  // PHObject virtual overloads
  
  virtual void         identify(std::ostream& os = std::cout) const {
    os << "---SvtxHit base class------------" << std::endl;
  }
  virtual void         Reset() {};
  virtual int          isValid() const {return 0;}
  virtual SvtxHit*     Clone() const {return NULL;}

  // digitized hit info
  
  virtual unsigned int get_id() const                        {return UINT_MAX;}
  virtual void         set_id(unsigned int id)               {}
  
  virtual unsigned int get_layer() const                     {return UINT_MAX;}
  virtual void         set_layer(unsigned int layer)         {}

  virtual unsigned int get_adc() const                       {return UINT_MAX;}
  virtual void         set_adc(unsigned int adc)             {}

  virtual float        get_e() const                         {return UINT_MAX;}
  virtual void         set_e(float e)                        {}
  
  virtual PHG4CellDefs::keytype get_cellid() const                    {return UINT_MAX;}
  virtual void         set_cellid(PHG4CellDefs::keytype cellid)       {}

protected:
  SvtxHit() {}
  
private:
  
  ClassDef(SvtxHit, 2);
};

#endif
