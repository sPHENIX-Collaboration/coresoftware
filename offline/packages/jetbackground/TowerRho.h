#ifndef JETBACKGROUND_TOWERRHO_H
#define JETBACKGROUND_TOWERRHO_H

//===========================================================
/// \file TowerRho.h
/// \brief PHObject to store rho and sigma for calorimeter towers on an event-by-event basis
/// \author Tanner Mengel
//===========================================================

#include <jetbase/Jet.h>

#include <phool/PHObject.h>

#include <iostream>

class TowerRho : public PHObject
{
 public:
  // enum for method of rho calculation
  enum Method
  {
    NONE = 0,
    AREA = 1,
    MULT = 2
  };

  ~TowerRho() override {};

  void identify(std::ostream &os = std::cout) const override { os << "TowerRho base class" << std::endl; };
  int isValid() const override { return 0; }

  // setters
  virtual void set_rho(float /*rho*/) {}
  virtual void set_sigma(float /*sigma*/) {}
  virtual void set_method(TowerRho::Method /*rho_method*/) {}

  // getters
  virtual float get_rho() { return 0; }
  virtual float get_sigma() { return 0; }
  virtual TowerRho::Method get_method() { return Method::NONE; }

 protected:
  TowerRho() {}  // ctor

 private:
  ClassDefOverride(TowerRho, 1);
};

#endif  // JETBACKGROUND_TOWERRHO_H
