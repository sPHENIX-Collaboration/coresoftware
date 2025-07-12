/*!
 * \file EventRho.h
 * \brief PHObject to store rho and sigma for calorimeter towers on an event-by-event basis. Options for rho calculation are AREA and MULT or NONE.
 * \author Tanner Mengel <tmengel@bnl.gov>
 * \version $Verison: 2.0.1 $
 * \date $Date: 02/01/2024. Revised 09/19/2024$
 */

#ifndef JETBACKGROUND_EVENTRHO_H
#define JETBACKGROUND_EVENTRHO_H

#include <jetbase/Jet.h>

#include <phool/PHObject.h>

#include <iostream>

class EventRho : public PHObject
{
 public:
  // enum for method of rho calculation
  enum Method
  {
    NONE = 0,
    AREA = 1,
    MULT = 2
  };

  ~EventRho() override {};

  void identify(std::ostream &os = std::cout) const override { os << "EventRho base class" << std::endl; };
  int isValid() const override { return 0; }

  // setters
  virtual void set_rho(float /*rho*/) {}
  virtual void set_sigma(float /*sigma*/) {}
  virtual void set_method(EventRho::Method /*rho_method*/) {}

  // getters
  virtual float get_rho() { return 0; }
  virtual float get_sigma() { return 0; }
  virtual EventRho::Method get_method() { return Method::NONE; }

 protected:
  EventRho() {}  // ctor

 private:
  ClassDefOverride(EventRho, 1);
};

#endif  // JETBACKGROUND_EVENTRHO_H
