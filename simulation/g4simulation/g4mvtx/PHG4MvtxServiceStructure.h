#ifndef G4MVTX_SERVICESTRUCTURE_H
#define G4MVTX_SERVICESTRUCTURE_H

#include <string>

class PHG4MvtxServiceStructure
{
 public:
  PHG4MvtxServiceStructure();

  explicit PHG4MvtxServiceStructure(const std::string &name,
                                    const float &thickness_copper,
                                    const float &thickness_carbon,
                                    const float &thickness_plastic,
                                    const float &zSouth,
                                    const float &zNorth,
                                    const float &rSouth,
                                    const float &rNorth);

  virtual ~PHG4MvtxServiceStructure(){};

  std::string get_name();
  float get_thickness_copper();
  float get_thickness_carbon();
  float get_thickness_plastic();
  float get_zSouth();
  float get_zNorth();
  float get_rSouth();
  float get_rNorth();

 private:
  const std::string m_name = "service";
  const float m_thickness_copper = 0.0;
  const float m_thickness_carbon = 0.0;
  const float m_thickness_plastic = 0.0;
  const float m_zSouth = 0.0;
  const float m_zNorth = 0.0;
  const float m_rSouth = 0.0;
  const float m_rNorth = 0.0;
};

#endif
