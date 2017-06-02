#ifndef __TOWERBACKGROUND_H__
#define __TOWERBACKGROUND_H__

#include <phool/PHObject.h>

#include <vector>

class TowerBackground : public PHObject
{
 public:
  virtual ~TowerBackground() {};

  virtual void identify(std::ostream &os = std::cout) const { os << "TowerBackground base class" << std::endl; };
  virtual void Reset() {}
  virtual int isValid() const { return 0; }

  virtual void set_UE( int layer, const std::vector<float> & UE ) {}
  virtual void set_v2( int layer, float v2 ) {}
  virtual void set_Psi2( int layer, float Psi2 ) {}
  
  virtual std::vector<float> get_UE( int layer ) { return std::vector<float>(); } ;
  virtual float get_v2( int layer ) { return 0; }
  virtual float get_Psi2( int layer ) { return 0; }

 protected:
  TowerBackground() {}

 private:

  ClassDef(TowerBackground, 1);
};

#endif  // __TOWERBACKGROUND_H__
