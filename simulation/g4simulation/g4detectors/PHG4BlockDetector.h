#ifndef PHG4BlockDetector_h
#define PHG4BlockDetector_h

#include "PHG4Parameters.h"

#include <g4main/PHG4Detector.h>

#include <Geant4/globals.hh>
#include <Geant4/G4Region.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Types.hh>

#include <map>


class G4Material;
class G4Box;
class G4LogicalVolume;
class G4Region;
class G4VPhysicalVolume;

class PHG4BlockDetector: public PHG4Detector
{

  public:

  //! constructor
  PHG4BlockDetector( PHCompositeNode *Node, PHG4Parameters *parameters, const std::string &dnam="BLOCK", const int lyr = 0 );

  //! destructor
  virtual ~PHG4BlockDetector( void )
  {}

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  //!@name volume accessors
  //@{
  bool IsInBlock(G4VPhysicalVolume*) const;
  //@}

  void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() const {return superdetector;}
  int get_Layer() const {return layer;}

  private:

  PHG4Parameters *params;
 
  G4VPhysicalVolume* block_physi;


  int layer;
  std::string superdetector;
  
};

#endif
