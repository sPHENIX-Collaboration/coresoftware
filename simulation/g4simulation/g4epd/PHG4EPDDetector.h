// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/* vim: set sw=2 ft=cpp: */

#ifndef G4EPD_PHG4EPDDETECTOR_H
#define G4EPD_PHG4EPDDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <cstdint>
#include <map>
#include <set>
#include <string>

class G4ExtrudedSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4EPDDisplayAction;
class PHG4Subsystem;
class PHParameters;

class PHG4EPDDetector : public PHG4Detector
{
 public:
  PHG4EPDDetector(PHG4Subsystem* subsys,
                  PHCompositeNode* node,
                  PHParameters* parameters,
                  std::string const& name);

  void ConstructMe(G4LogicalVolume* world) override;

  int IsInDetector(G4VPhysicalVolume*) const;

  uint32_t module_id_for(uint32_t index, uint32_t slice, uint32_t side);
  uint32_t module_id_for(G4VPhysicalVolume* volume);

  void SuperDetector(std::string const& name) { superdetector = name; }
  const std::string SuperDetector() const { return superdetector; }

  PHG4EPDDisplayAction* GetDisplayAction() { return m_DisplayAction; }

 private:
  G4ExtrudedSolid* construct_block(int32_t index);

  PHG4EPDDisplayAction* m_DisplayAction = nullptr;
  PHParameters* m_Params = nullptr;

  int m_ActiveFlag = 0;
  int m_SupportActiveFlag = 0;

  std::set<G4LogicalVolume*> m_SupportLogVolSet;
  std::set<G4LogicalVolume*> m_ActiveLogVolSet;

  std::map<G4VPhysicalVolume*, uint32_t> m_volumes;

  std::string superdetector;
};

#endif /* G4EPD_PHG4EPDETECTOR_H */
