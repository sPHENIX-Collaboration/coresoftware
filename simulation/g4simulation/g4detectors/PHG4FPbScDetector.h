// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FPBSCDETECTOR_H
#define G4DETECTORS_PHG4FPBSCDETECTOR_H

#include "g4main/PHG4Detector.h"

#include <Geant4/G4Region.hh>
#include <Geant4/G4String.hh>           // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Types.hh>

#include <cstddef>                     // for size_t
#include <map>
#include <string>                       // for string


class G4Material;
class G4Box;
class G4LogicalVolume;
class G4UserSteppingAction;
class G4VPhysicalVolume;
class PHCompositeNode;

class PHG4FPbScDetector: public PHG4Detector
{
  public:
    
  PHG4FPbScDetector( PHCompositeNode *Node, const std::string &nam );
    
    virtual ~PHG4FPbScDetector( void ){}
    
    virtual void Construct( G4LogicalVolume* world );
    
    virtual G4UserSteppingAction* GetSteppingAction() 
    { 
      if ( _region )
        return _region->GetRegionalSteppingAction();
      else return 0;
    }
    
    bool isInScintillator(G4VPhysicalVolume * volume);
    int getScintillatorLayer(G4VPhysicalVolume * volume);
    // compute tower index
    unsigned int computeIndex(unsigned int layer, G4double x, G4double y, G4double z, G4double& xcenter, G4double& ycenter, G4double& zcenter);
    void set_Place(G4double x, G4double y, G4double z)
    {
      x_position = x * cm;
      y_position = y * cm;
      z_position = z * cm;
    }

  private:
    
    G4double tower_cross_section;
    size_t segments_per_column, segments_per_height;
    G4double length, height, absorber_thickness, scintillator_thickness;
    unsigned int nlayers;//, segments_per_thickness;
    G4double x_position, y_position, z_position;
    G4double layer_separation;
    
    G4Material* SetMaterial(G4String);
    
    G4Material* AbsorberMaterial;
    G4Material* ScintillatorMaterial;
    
    std::map<unsigned int, G4Box* > absorber_solid_;
    std::map<unsigned int, G4LogicalVolume* > absorber_logic_;
    std::map<unsigned int, G4VPhysicalVolume* > absorber_physi_;
    
    std::map<unsigned int, G4Box* > scintillator_solid_;
    std::map<unsigned int, G4LogicalVolume* > scintillator_logic_;
    std::map<unsigned int, G4VPhysicalVolume* > scintillator_physi_;
    
    G4Region* _region;
};


#endif
