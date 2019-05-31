// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BLOCKDISPLAYACTION_H
#define G4DETECTORS_PHG4BLOCKDISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

#include <string>                      // for string

class G4VisAttributes;
class G4LogicalVolume;
class G4VPhysicalVolume;
class PHParameters;

class PHG4BlockDisplayAction : public PHG4DisplayAction
{
 public:
  PHG4BlockDisplayAction(const std::string &name, PHParameters *parameters);

  virtual ~PHG4BlockDisplayAction();

  void ApplyDisplayAction(G4VPhysicalVolume *physvol);
  void SetMyVolume(G4LogicalVolume *vol) { m_MyVolume = vol; }

 private:
  PHParameters *m_Params;
  G4LogicalVolume *m_MyVolume;
  G4VisAttributes *m_VisAtt;
};

#endif  // G4DETECTORS_PHG4BLOCKDISPLAYACTION_H
