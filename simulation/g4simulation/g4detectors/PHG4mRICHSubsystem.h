/*===============================================================*
 *                        March 2nd 2017                         *
 *        mRICH Subsystem created by Cheuk-Ping Wong @GSU        *
 *===============================================================*/
#ifndef G4DETECTORS_PHG4MRICHSUBSYSTEM_H
#define G4DETECTORS_PHG4MRICHSUBSYSTEM_H

#include "PHG4DetectorSubsystem.h"

class PHG4mRICHDetector;
class PHG4BlockSteppingAction;
class PHG4EventAction;

class PHG4mRICHSubsystem: public PHG4DetectorSubsystem
{
 public:

  //! constructor
  PHG4mRICHSubsystem( const std::string &name = "BLOCK", const int layer = 0 );

  //! destructor
  virtual ~PHG4mRICHSubsystem( void )
  {}

  int InitSubsystem(PHCompositeNode *);

  //! InitRunSubsystem
  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode *);

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode *);

  //! accessors (reimplemented)
  virtual PHG4Detector* GetDetector( void ) const;
  virtual PHG4SteppingAction* GetSteppingAction( void ) const {return  _steppingAction;}

  PHG4EventAction* GetEventAction() const {return _eventAction;}

 private:
  void SetDefaultParameters();       //set external parameter

  //! detector geometry
  /*! defives from PHG4Detector */
  //PHG4BlockDetector* _detector;
  //int _single_mRICH;
  PHG4mRICHDetector* _detector;
  std::string _detectorName;
//  int layer;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* _steppingAction;
  PHG4EventAction *_eventAction;

};

#endif
