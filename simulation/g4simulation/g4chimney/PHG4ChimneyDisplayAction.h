#ifndef PHG4CHIMNEYDISPLAYACTION_H
#define PHG4CHIMNEYDISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

#include <map>
#include <string>
#include <vector>

class G4LogicalVolume;
class G4VisAttributes;
class G4VPhysicalVolume;

class PHG4ChimneyDisplayAction : public PHG4DisplayAction
{
 public:
  explicit PHG4ChimneyDisplayAction(const std::string &name);
  ~PHG4ChimneyDisplayAction() override;
  void ApplyDisplayAction(G4VPhysicalVolume *physvol) override;
  void AddVolume(G4LogicalVolume *logvol, const std::string &mat) { m_LogicalVolumeMap[logvol] = mat; }

 private:
  std::map<G4LogicalVolume *, std::string> m_LogicalVolumeMap;
  std::vector<G4VisAttributes *> m_VisAttVec;
};

#endif
