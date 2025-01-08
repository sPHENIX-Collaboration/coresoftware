#ifndef G4MVTX_CABLE_H
#define G4MVTX_CABLE_H

#include <string>

class PHG4MvtxCable
{
 public:
  PHG4MvtxCable();

  explicit PHG4MvtxCable(const std::string &name,
                         const std::string &coreMaterial,
                         const double &coreRadius,
                         const double &sheathRadius,
                         const double &xSouth, const double &xNorth,
                         const double &ySouth, const double &yNorth,
                         const double &zSouth, const double &zNorth,
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

  virtual ~PHG4MvtxCable() = default;

  std::string get_name() { return m_name; };
  std::string get_coreMaterial() { return m_coreMaterial; };
  double get_coreRadius() { return m_coreRadius; };
  double get_sheathRadius() { return m_sheathRadius; };
  double get_xSouth() { return m_xSouth; };
  double get_xNorth() { return m_xNorth; };
  double get_ySouth() { return m_ySouth; };
  double get_yNorth() { return m_yNorth; };
  double get_zSouth() { return m_zSouth; };
  double get_zNorth() { return m_zNorth; };
  std::string get_color() { return m_color; };

 private:
  const std::string m_name = "cable";
  const std::string m_coreMaterial = "G4_Cu";
  const double m_coreRadius = 1;
  const double m_sheathRadius = 2;
  const double m_xSouth = 0.;
  const double m_xNorth = 1.;
  const double m_ySouth = 0.;
  const double m_yNorth = 1.;
  const double m_zSouth = 0.;
  const double m_zNorth = 1.;
  const std::string m_color = "red";
};

#endif
