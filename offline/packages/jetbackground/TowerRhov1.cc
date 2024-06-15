#include "TowerRhov1.h"

#include <ostream>
#include <cstdlib> // for exit


TowerRhov1::TowerRhov1()
  : m_tower_rho(0)
  , m_tower_sigma(0)
  , m_rho_method_type(TowerRho::Method::NONE)
{
  
}

void TowerRhov1::identify(std::ostream& os) const
{
  os << "TowerRhov1: " << std::endl;
  os << "\t method: " << get_method_string(m_rho_method_type) << std::endl;
  os << "\trho = " << m_tower_rho << ", sigma(rho) = " << m_tower_sigma << std::endl;
  os << "===============================";  
  return ;
}

void TowerRhov1::set_method(TowerRho::Method method)
{
  get_method_string(method); // check if method is valid
  m_rho_method_type = method;
  return ;
}

std::string TowerRhov1::get_method_string(TowerRho::Method method)
{
  switch (method)
  {
  case TowerRho::Method::NONE:
    return "NONE";
    break;
  case TowerRho::Method::AREA:
    return "AREA";
    break;
  case TowerRho::Method::MULT:
    return "MULT";
    break;
  default:
    std::cout << "ERROR: rho method not recognized" << std::endl;
    std::cout << "rho method must be 1 (area) or 2 (mult)" << std::endl;
    exit(-1);
  }
  return "NONE";
}


