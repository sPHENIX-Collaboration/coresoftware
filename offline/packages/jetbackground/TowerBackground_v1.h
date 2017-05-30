#ifndef __TOWERBACKGROUND_V1_H__
#define __TOWERBACKGROUND_V1_H__

#include "TowerBackground.h"

#include <phool/PHObject.h>

class TowerBackground_v1 : public TowerBackground
{
 public:
  TowerBackground_v1();
  virtual ~TowerBackground_v1();

  void identify(std::ostream &os = std::cout) const;
  void Reset() {}
  int isValid() const { return 1; }

  virtual void set_UE( int layer, const std::vector<float> & UE ) { _UE[layer] = UE; }
  virtual void set_v2( int layer, float v2 ) { _v2[layer] = v2; }
  virtual void set_Psi2( int layer, float Psi2 ) { _Psi2[layer] = Psi2; }
  
  virtual std::vector<float> get_UE( int layer ) { return _UE[layer]; }
  virtual float get_v2( int layer ) { return _v2[layer]; }
  virtual float get_Psi2( int layer ) { return _Psi2[layer]; }

 private:

  std::vector< std::vector< float> > _UE;
  float _v2[3];
  float _Psi2[3];

  ClassDef(TowerBackground_v1, 1);
};

#endif  // __TOWERBACKGROUND_V1_H__
