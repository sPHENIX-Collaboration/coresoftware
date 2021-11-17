// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CONEDISPLAYACTION_H
#define G4DETECTORS_PHG4CONEDISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

#include <string>  // for string

class G4Colour;
class G4LogicalVolume;
class G4VisAttributes;
class G4VPhysicalVolume;
class PHParameters;

class PHG4ConeDisplayAction : public PHG4DisplayAction
{
 public:
  PHG4ConeDisplayAction(const std::string &name, PHParameters *parameters);

  ~PHG4ConeDisplayAction() override;

  void ApplyDisplayAction(G4VPhysicalVolume *physvol) override;
  void SetMyVolume(G4LogicalVolume *vol) { m_MyVolume = vol; }
  void SetColor(const double red, const double green, const double blue, const double alpha = 1.);

 private:
  PHParameters *m_Params;
  G4LogicalVolume *m_MyVolume;
  G4VisAttributes *m_VisAtt;
  G4Colour *m_Colour;
};

#endif  // G4DETECTORS_PHG4CONEDISPLAYACTION_H
