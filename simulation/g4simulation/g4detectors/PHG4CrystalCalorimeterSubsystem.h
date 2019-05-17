// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CRYSTALCALORIMETERSUBSYSTEM_H
#define G4DETECTORS_PHG4CRYSTALCALORIMETERSUBSYSTEM_H

#include <g4main/PHG4Subsystem.h>

class PHG4CrystalCalorimeterDetector;
class PHG4DisplayAction;
class PHG4ProjCrystalCalorimeterDetector;
class PHG4SteppingAction;

class PHG4CrystalCalorimeterSubsystem: public PHG4Subsystem
{

public:

  /** Constructor
   */
  PHG4CrystalCalorimeterSubsystem( const std::string &name = "CRYSTAL_DEFAULT", const int layer = 0 );

  /** Destructor
   */
  virtual ~PHG4CrystalCalorimeterSubsystem( );

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
  PHG4Detector* GetDetector( void ) const;
  PHG4SteppingAction* GetSteppingAction(  ) const { return m_SteppingAction; }
  PHG4DisplayAction *GetDisplayAction() const { return m_DisplayAction; }

  /** Set mapping file for calorimeter towers
   */
  void SetTowerMappingFile( const std::string &filename )
  {
    mappingfile_ = filename;
  }

  /** Select projective geometry for calorimeter
   */
  void SetProjectiveGeometry( const std::string &filename1 , const std::string &filename2 ) {
    mappingfile_ = filename1;
    mappingfile_4x4_construct_ = filename2;
    projective_ = true;
  }
//  void SetGeometryConfiguration

  /** Enum for different geometry configurations
   */
//  static const enum GeometryConfiguration {
//    kNonProjective,
//    kProjectiveV1
//  };

private:


  //  GeometryConfiguration current_geom_config_;

  /** Pointer to the Geant4 implementation of the detector
   */
  PHG4CrystalCalorimeterDetector* detector_;

  /** Stepping action
   */
  PHG4SteppingAction *m_SteppingAction;
  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction *m_DisplayAction;

  int active;

  std::string detector_type;
  std::string mappingfile_;
  std::string mappingfile_4x4_construct_;

  bool projective_;

};

#endif
