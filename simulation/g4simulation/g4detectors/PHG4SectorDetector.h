#ifndef PHG4SectorDetector_h
#define PHG4SectorDetector_h

#include "g4main/PHG4Detector.h"

#include <Geant4/G4Types.hh>
#include <Geant4/G4Region.hh>
#include <Geant4/globals.hh>

#include "PHG4SectorConstructor.h"

#include <map>


class G4Material;
class G4Cons;
class G4LogicalVolume;
class G4Region;
class G4VPhysicalVolume;

class PHG4SectorDetector: public PHG4Detector, public PHG4Sector::PHG4SectorConstructor
{

  public:

  //! constructor
  PHG4SectorDetector( PHCompositeNode *Node, const std::string &dnam="SECTOR");

  //! destructor
  virtual ~PHG4SectorDetector( void )
  {}

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  //!@name volume accessors
  //@{
  bool IsInSectorActive(G4VPhysicalVolume*);
  bool IsInSectorInactive(G4VPhysicalVolume*);
  //@}


  virtual G4UserSteppingAction* GetSteppingAction() 
  { 
    if ( _region )
      return _region->GetRegionalSteppingAction();
    else return 0;
  }


  void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() const {return superdetector;}


  virtual void OverlapCheck(const bool chk = true)
  {
    PHG4Detector::OverlapCheck(chk);
    PHG4SectorConstructor::OverlapCheck(chk);
  }

  private:


  G4Region* _region;

  std::string superdetector;
  
};

#endif
