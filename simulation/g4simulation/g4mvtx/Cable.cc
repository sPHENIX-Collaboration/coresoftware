#include "Cable.h"

Cable::Cable(const std::string &name,
             const std::string &coreMaterial,
             const float &coreRadius,
             const float &sheathRadius,
             const float &xSouth, const float &xNorth,
             const float &ySouth, const float &yNorth,
             const float &zSouth, const float &zNorth,
             const std::vector<float> &RGB)
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
  , m_RGB(RGB)
{
}

std::string Cable::get_name() { return m_name; }
std::string Cable::get_coreMaterial() { return m_coreMaterial; }
float Cable::get_coreRadius() { return m_coreRadius; }
float Cable::get_sheathRadius() { return m_sheathRadius; }
float Cable::get_xSouth() { return m_xSouth; }
float Cable::get_xNorth() { return m_xNorth; }
float Cable::get_ySouth() { return m_ySouth; }
float Cable::get_yNorth() { return m_yNorth; }
float Cable::get_zSouth() { return m_zSouth; }
float Cable::get_zNorth() { return m_zNorth; }
std::vector<float> Cable::get_RGB() { return m_RGB; }
