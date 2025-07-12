#include "EventRhov1.h"

#include <cstdlib>

EventRhov1::EventRhov1()
  : m_tower_rho(0)
  , m_tower_sigma(0)
  , m_rho_method_type(EventRho::Method::NONE)
{
}

void EventRhov1::identify(std::ostream& os) const
{
  os << "EventRhov1: " << std::endl;
  os << "\t method: " << get_method_string(m_rho_method_type) << std::endl;
  os << "\trho = " << m_tower_rho << ", sigma(rho) = " << m_tower_sigma << std::endl;
  os << "===============================";
  return;
}

void EventRhov1::set_method(EventRho::Method method)
{
  get_method_string(method);
  m_rho_method_type = method;
  return;
}

std::string EventRhov1::get_method_string(EventRho::Method method)
{
  switch (method)
  {
  case EventRho::Method::NONE:
    return "NONE";
    break;
  case EventRho::Method::AREA:
    return "AREA";
    break;
  case EventRho::Method::MULT:
    return "MULT";
    break;
  default:
    std::cout << "ERROR: rho method not recognized" << std::endl;
    std::cout << "rho method must be 1 (area) or 2 (mult)" << std::endl;
    exit(-1);
  }
  return "NONE";
}
