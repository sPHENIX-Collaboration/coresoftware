#include "ServiceStructure.h"

ServiceStructure::ServiceStructure(const std::string &name,
                                   const float &thickness_copper,
                                   const float &thickness_carbon,
                                   const float &zSouth,
                                   const float &zNorth,
                                   const float &rSouth,
                                   const float &rNorth)
  : m_name(name)
  , m_thickness_copper(thickness_copper)
  , m_thickness_carbon(thickness_carbon)
  , m_zSouth(zSouth)
  , m_zNorth(zNorth)
  , m_rSouth(rSouth)
  , m_rNorth(rNorth)
{
}

std::string ServiceStructure::get_name() { return m_name; }
float ServiceStructure::get_thickness_copper() { return m_thickness_copper; }
float ServiceStructure::get_thickness_carbon() { return m_thickness_carbon; }
float ServiceStructure::get_zSouth() { return m_zSouth; }
float ServiceStructure::get_zNorth() { return m_zNorth; }
float ServiceStructure::get_rSouth() { return m_rSouth; }
float ServiceStructure::get_rNorth() { return m_rNorth; }
