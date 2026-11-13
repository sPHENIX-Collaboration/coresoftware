// sPHENIX magnet chimney (valve box extension): box on the cryostat, horizontal Al tube, copper bus. Draft 2026-10-26.
#ifndef PHG4CHIMNEYDETECTOR_H
#define PHG4CHIMNEYDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <set>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4ChimneyDisplayAction;
class PHG4Subsystem;
class PHParameters;

class PHG4ChimneyDetector : public PHG4Detector
{
 public:
  PHG4ChimneyDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);
  ~PHG4ChimneyDetector() override = default;

  void ConstructMe(G4LogicalVolume *world) override;
  int IsInChimney(G4VPhysicalVolume *) const;

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string &SuperDetector() const { return m_SuperDetector; }
  PHG4ChimneyDisplayAction *GetDisplayAction() { return m_DisplayAction; }

 private:
  PHG4ChimneyDisplayAction *m_DisplayAction = nullptr;
  PHParameters *m_Params = nullptr;
  std::string m_SuperDetector;
  std::set<G4LogicalVolume *> m_LogicalVolSet;
};

#endif
