#ifndef G4MVTX_CABLE_H
#define G4MVTX_CABLE_H

#include <string>

class PHG4MvtxCable
{
 public:
  PHG4MvtxCable();

  explicit PHG4MvtxCable(const std::string &name,
                         const std::string &coreMaterial,
                         const float &coreRadius,
                         const float &sheathRadius,
                         const float &xSouth, const float &xNorth,
                         const float &ySouth, const float &yNorth,
                         const float &zSouth, const float &zNorth,
                         const std::string &color);

  virtual ~PHG4MvtxCable(){};

  std::string get_name();
  std::string get_coreMaterial();
  float get_coreRadius();
  float get_sheathRadius();
  float get_xSouth();
  float get_xNorth();
  float get_ySouth();
  float get_yNorth();
  float get_zSouth();
  float get_zNorth();
  std::string get_color();

 private:
  const std::string m_name = "cable";
  const std::string m_coreMaterial = "G4_Cu";
  const float m_coreRadius = 1;
  const float m_sheathRadius = 2;
  const float m_xSouth = 0.;
  const float m_xNorth = 1.;
  const float m_ySouth = 0.;
  const float m_yNorth = 1.;
  const float m_zSouth = 0.;
  const float m_zNorth = 1.;
  const std::string m_color = "red";
};

#endif
