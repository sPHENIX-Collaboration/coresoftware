// $Id: $

/*!
 * \file PHG4GDMLDetector.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHG4GDMLDetector_H_
#define PHG4GDMLDetector_H_

#include <g4main/PHG4Detector.h>
#include <string>

#include <Geant4/G4Types.hh>
#include <Geant4/globals.hh>

class PHCompositeNode;
class PHParameters;
class G4UserSteppingAction;
class G4UserSteppingAction;

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
