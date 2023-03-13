// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MICROMEGAS_PHG4MICROMEGASDISPLAYACTION_H
#define G4MICROMEGAS_PHG4MICROMEGASDISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

#include <Geant4/G4Color.hh>

#include <map>
#include <string>  // for string
#include <vector>

class G4LogicalVolume;
class G4VisAttributes;
class G4VPhysicalVolume;
class G4Colour;

class PHG4MicromegasDisplayAction : public PHG4DisplayAction
{
 public:
  PHG4MicromegasDisplayAction(const std::string &name);

  ~PHG4MicromegasDisplayAction() override;

  void ApplyDisplayAction(G4VPhysicalVolume *physvol) override;
  void AddVolume(G4LogicalVolume *logvol, const G4Colour &col ) { m_LogicalVolumeMap[logvol] = col; }

 private:
  std::map<G4LogicalVolume *, G4Colour> m_LogicalVolumeMap;
  std::vector<G4VisAttributes *> m_VisAttVec;
};

#endif  // G4MICROMEGAS_PHG4MICROMEGASDISPLAYACTION_H
