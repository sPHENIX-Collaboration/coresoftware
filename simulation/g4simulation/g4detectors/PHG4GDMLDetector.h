// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $Id: $

/*!
 * \file PHG4GDMLDetector.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4DETECTORS_PHG4GDMLDETECTOR_H
#define G4DETECTORS_PHG4GDMLDETECTOR_H

#include <g4main/PHG4Detector.h>


#include <Geant4/G4Types.hh>

#include <string>

class G4AssemblyVolume;
class G4LogicalVolume;
class G4UserSteppingAction;
class PHCompositeNode;
class PHParameters;

/*!
 * \brief PHG4GDMLDetector is a generic detector built from a GDML import
 */
class PHG4GDMLDetector : public PHG4Detector
{
 public:
  PHG4GDMLDetector(PHCompositeNode* Node, const std::string& dnam, PHParameters* parameters);

  virtual ~PHG4GDMLDetector();

  //! construct
  void Construct(G4LogicalVolume* world);

  G4UserSteppingAction* GetSteppingAction()
  {
    return nullptr;
  }

  void Print(const std::string& what = "ALL") const;

 private:
  void SetDisplayProperty( G4AssemblyVolume* av);
  void SetDisplayProperty( G4LogicalVolume* lv);

  std::string m_GDMPath;
  std::string m_TopVolName;

  G4double m_placeX;
  G4double m_placeY;
  G4double m_placeZ;

  G4double m_rotationX;
  G4double m_rotationY;
  G4double m_rotationZ;
};

#endif /* PHG4GDMLDetector_H_ */
