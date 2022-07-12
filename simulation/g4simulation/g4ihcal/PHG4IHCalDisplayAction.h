// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4IHCAL_PHG4INNERHCALDISPLAYACTION_H
#define G4IHCAL_PHG4INNERHCALDISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

#include <set>
#include <string>  // for string
#include <vector>

class G4LogicalVolume;
class G4VisAttributes;
class G4VPhysicalVolume;

class PHG4IHCalDisplayAction : public PHG4DisplayAction
{
 public:
  explicit PHG4IHCalDisplayAction(const std::string &name);

  ~PHG4IHCalDisplayAction() override;

  void ApplyDisplayAction(G4VPhysicalVolume *physvol) override;
  void SetMyTopVolume(G4VPhysicalVolume *vol) { m_MyTopVolume = vol; }
  void AddScintiVolume(G4LogicalVolume *vol) { m_ScintiLogVolSet.insert(vol); }
  void AddSteelVolume(G4LogicalVolume *vol) { m_SteelVol = vol; }

 private:
  G4VPhysicalVolume *m_MyTopVolume;
  G4LogicalVolume *m_SteelVol;
  std::vector<G4VisAttributes *> m_VisAttVec;
  std::set<G4LogicalVolume *> m_ScintiLogVolSet;
};

#endif  // G4IHCAL_PHG4INNERHCALDISPLAYACTION_H
