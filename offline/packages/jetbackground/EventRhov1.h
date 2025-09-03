/*!
 * \file EventRhov1.h
 * \brief PHObject to store rho and sigma for calorimeter towers on an event-by-event basis. Options for rho calculation are AREA and MULT or NONE. Copy of EventRhov1 with a differnet name
 * \author Tanner Mengel <tmengel@bnl.gov>
 * \version $Verison: 2.0.1 $
 * \date $Date: 02/01/2024. $
 */

#ifndef JETBACKGROUND_EVENTRHOV1_H
#define JETBACKGROUND_EVENTRHOV1_H

#include "EventRho.h"

class EventRhov1 : public EventRho
{
 public:
  EventRhov1();
  ~EventRhov1() override = default;

  void identify(std::ostream &os = std::cout) const override;
  void Reset() override {}  // new in v1
  int isValid() const override { return 1; }

  // setters
  void set_rho(float rho) override { m_tower_rho = rho; }
  void set_sigma(float sigma) override { m_tower_sigma = sigma; }
  void set_method(EventRho::Method method) override;

  // getters
  float get_rho() override { return m_tower_rho; }
  float get_sigma() override { return m_tower_sigma; }
  EventRho::Method get_method() override { return m_rho_method_type; }

  // static method to string conversion
  static std::string get_method_string(EventRho::Method method);  // new in v1

 private:
  float m_tower_rho{0};                                        // momentum density
  float m_tower_sigma{0};                                      // sigma of momentum density
  EventRho::Method m_rho_method_type{EventRho::Method::NONE};  // method of rho calculation

  ClassDefOverride(EventRhov1, 1);
};

#endif  // JETBACKGROUND_EVENTRHOV1_H
