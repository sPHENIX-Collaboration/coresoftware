// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FORWARDECALSUBSYSTEM_H
#define G4DETECTORS_PHG4FORWARDECALSUBSYSTEM_H

#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

class PHG4DisplayAction;
class PHG4ForwardEcalDetector;
class PHG4ForwardEcalSteppingAction;

class PHG4ForwardEcalSubsystem: public PHG4Subsystem
{

public:

  /** Constructor
   */
  PHG4ForwardEcalSubsystem( const std::string &name = "FORWARD_ECAL_DEFAULT", const int layer = 0 );

  /** Destructor
   */
  virtual ~PHG4ForwardEcalSubsystem();


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
  PHG4SteppingAction* GetSteppingAction( void ) const;

   PHG4DisplayAction* GetDisplayAction() const { return m_DisplayAction; }

 /** Set mapping file for calorimeter towers
   */
  void SetTowerMappingFile( const std::string &filename )
  {
    mappingfile_ = filename;
  }

  void SetActive(const int i = 1){active = i;}
  void SetAbsorberActive(const int i = 1){absorber_active = i;}
  void BlackHole(const int i=1){blackhole = i;}

  void SetEICDetector(){EICDetector = 1;}
  void SetfsPHENIXDetector(){EICDetector = 0;}

private:

  /** Pointer to the Geant4 implementation of the detector
   */
  PHG4ForwardEcalDetector* detector_;

  /** Stepping action
   */
  PHG4ForwardEcalSteppingAction* steppingAction_;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction;

  int active;
  int absorber_active; 
  int blackhole; 

  std::string detector_type;
  std::string mappingfile_;

  int EICDetector; 

};

#endif
