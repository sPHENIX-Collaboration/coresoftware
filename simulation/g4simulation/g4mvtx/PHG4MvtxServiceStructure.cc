#include "PHG4MvtxServiceStructure.h"

PHG4MvtxServiceStructure::PHG4MvtxServiceStructure(const std::string &name,
                                                   const float &thickness_copper,
                                                   const float &thickness_carbon,
                                                   const float &thickness_plastic,
                                                   const float &zSouth,
                                                   const float &zNorth,
                                                   const float &rSouth,
                                                   const float &rNorth)
  : m_name(name)
  , m_thickness_copper(thickness_copper)
  , m_thickness_carbon(thickness_carbon)
  , m_thickness_plastic(thickness_plastic)
  , m_zSouth(zSouth)
  , m_zNorth(zNorth)
  , m_rSouth(rSouth)
  , m_rNorth(rNorth)
{
}

std::string PHG4MvtxServiceStructure::get_name() { return m_name; }
float PHG4MvtxServiceStructure::get_thickness_copper() { return m_thickness_copper; }
float PHG4MvtxServiceStructure::get_thickness_carbon() { return m_thickness_carbon; }
float PHG4MvtxServiceStructure::get_thickness_plastic() { return m_thickness_plastic; }
float PHG4MvtxServiceStructure::get_zSouth() { return m_zSouth; }
float PHG4MvtxServiceStructure::get_zNorth() { return m_zNorth; }
float PHG4MvtxServiceStructure::get_rSouth() { return m_rSouth; }
float PHG4MvtxServiceStructure::get_rNorth() { return m_rNorth; }
