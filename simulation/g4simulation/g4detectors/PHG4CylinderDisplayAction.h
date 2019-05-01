// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4MAIN_PHG4CYLINDERDISPLAYACTION_H
#define G4MAIN_PHG4CYLINDERDISPLAYACTION_H

#include <g4main/PHG4DisplayAction.h>

class G4VPhysicalVolume;
class PHParameters;

class PHG4CylinderDisplayAction: public PHG4DisplayAction
{
public:
  PHG4CylinderDisplayAction( const std::string &name,  PHParameters *parameters );

  virtual ~PHG4CylinderDisplayAction() {}

  void ApplyDisplayAction(G4VPhysicalVolume* physvol);
  void SetVolName(const std::string &name) {volname = name;}
  std::string GetVolName() const {return volname;}
  void SetMyVolume(G4VPhysicalVolume *vol) {m_MyVolume = vol;}

protected:
  bool CheckVolume(G4VPhysicalVolume *physvol);
  void ApplyVisAttributes(G4VPhysicalVolume *vol);

private:
  std::string volname;

  PHParameters *m_Params;
  G4VPhysicalVolume *m_MyVolume;
};


#endif // G4MAIN_PHG4CYLINDERDISPLAYACTION_H
