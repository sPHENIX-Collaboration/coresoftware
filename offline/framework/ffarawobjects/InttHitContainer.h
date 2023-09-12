#ifndef FUN4ALLRAW_INTTHITCONTAINER_H
#define FUN4ALLRAW_INTTHITCONTAINER_H

#include <phool/PHObject.h>

class InttHit;

class  InttHitContainer: public PHObject  
{
public:
  InttHitContainer() = default;
  virtual ~InttHitContainer() = default;

  virtual InttHit *AddHit() {return nullptr;}
  virtual unsigned int get_nhits() {return 0;}
  virtual InttHit *get_hit(unsigned int) {return nullptr;}

private:
  ClassDefOverride(InttHitContainer,1)    
};

#endif
