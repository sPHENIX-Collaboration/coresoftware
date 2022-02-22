#ifndef G4MVTX_SERVICESTRUCTURE_H
#define G4MVTX_SERVICESTRUCTURE_H

#include <string>

class ServiceStructure
{
 public:
  ServiceStructure();

  explicit ServiceStructure(const std::string &name,
                            const float &thickness_copper,
                            const float &thickness_carbon,
                            const float &zSouth,
                            const float &zNorth,
                            const float &rSouth,
                            const float &rNorth);

  virtual ~ServiceStructure(){};

  std::string get_name();
  float get_thickness_copper();
  float get_thickness_carbon();
  float get_zSouth();
  float get_zNorth();
  float get_rSouth();
  float get_rNorth();

 private:
  const std::string m_name = "service";
  const float m_thickness_copper = 0.0;
  const float m_thickness_carbon = 0.0;
  const float m_zSouth = 0.0;
  const float m_zNorth = 0.0;
  const float m_rSouth = 0.0;
  const float m_rNorth = 0.0;
};

#endif
