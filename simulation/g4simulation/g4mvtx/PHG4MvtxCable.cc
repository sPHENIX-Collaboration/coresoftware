#include "PHG4MvtxCable.h"

PHG4MvtxCable::PHG4MvtxCable(const std::string &name,
                             const std::string &coreMaterial,
                             const float &coreRadius,
                             const float &sheathRadius,
                             const float &xSouth, const float &xNorth,
                             const float &ySouth, const float &yNorth,
                             const float &zSouth, const float &zNorth,
                             const std::string &color)
  : m_name(name)
  , m_coreMaterial(coreMaterial)
  , m_coreRadius(coreRadius)
  , m_sheathRadius(sheathRadius)
  , m_xSouth(xSouth)
  , m_xNorth(xNorth)
  , m_ySouth(ySouth)
  , m_yNorth(yNorth)
  , m_zSouth(zSouth)
  , m_zNorth(zNorth)
  , m_color(color)
{
}

std::string PHG4MvtxCable::get_name() { return m_name; }
std::string PHG4MvtxCable::get_coreMaterial() { return m_coreMaterial; }
float PHG4MvtxCable::get_coreRadius() { return m_coreRadius; }
float PHG4MvtxCable::get_sheathRadius() { return m_sheathRadius; }
float PHG4MvtxCable::get_xSouth() { return m_xSouth; }
float PHG4MvtxCable::get_xNorth() { return m_xNorth; }
float PHG4MvtxCable::get_ySouth() { return m_ySouth; }
float PHG4MvtxCable::get_yNorth() { return m_yNorth; }
float PHG4MvtxCable::get_zSouth() { return m_zSouth; }
float PHG4MvtxCable::get_zNorth() { return m_zNorth; }
std::string PHG4MvtxCable::get_color() { return m_color; }
