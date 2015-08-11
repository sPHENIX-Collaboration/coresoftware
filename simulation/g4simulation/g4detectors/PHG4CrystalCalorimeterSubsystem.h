#ifndef PHG4CrystalCalorimeterSubsystem_h
#define PHG4CrystalCalorimeterSubsystem_h

#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

class PHG4CrystalCalorimeterDetector;
class PHG4CrystalCalorimeterSteppingAction;
class PHG4EventAction;

class PHG4CrystalCalorimeterSubsystem: public PHG4Subsystem
{

public:

  /** Constructor
   */
  PHG4CrystalCalorimeterSubsystem( const std::string &name = "CRYSTAL_DEFAULT", const int layer = 0 );

  /** Destructor
   */
  virtual ~PHG4CrystalCalorimeterSubsystem( void )
  {}

  /**
     Creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
     Creates the stepping action and place it on the node tree, under "ACTIONS" node
     Creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int Init(PHCompositeNode *);

  /** Event processing
   */
  int process_event(PHCompositeNode *);

  /** Accessors (reimplemented)
   */
  virtual PHG4Detector* GetDetector( void ) const;
  virtual PHG4SteppingAction* GetSteppingAction( void ) const;

  /** Set mapping file for calorimeter towers
   */
  void SetTowerMappingFile( std::string filename )
  {
    mappingfile_ = filename;
  }

  /** Select mapping file for calorimeter tower
   */
  void SetProjectiveGeometry( std::string filename1 , std::string filename2 ) {
    mappingfile_ = filename1;
    mappingfile_4x4_construct_ = filename2;
    projective_ = true;
  }


private:

  /** Pointer to the Geant4 implementation of the detector
   */
  PHG4CrystalCalorimeterDetector* detector_;

  /** Stepping action
   */
  PHG4CrystalCalorimeterSteppingAction* steppingAction_;
  PHG4EventAction *eventAction_;

  G4String material;
  int active;

  std::string detector_type;
  std::string mappingfile_;
  std::string mappingfile_4x4_construct_;
  G4bool projective_;

};

#endif
