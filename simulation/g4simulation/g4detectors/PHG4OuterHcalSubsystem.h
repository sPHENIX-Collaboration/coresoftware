#ifndef PHG4OuterHcalSubsystem_h
#define PHG4OuterHcalSubsystem_h

#include "PHG4DetectorSubsystem.h"

#include <map>
#include <string>

class PHG4OuterHcalDetector;
class PHG4OuterHcalSteppingAction;

class PHG4OuterHcalSubsystem: public PHG4DetectorSubsystem
{

  public:

  //! constructor
  PHG4OuterHcalSubsystem( const std::string &name = "HCALOUT", const int layer = 0 );

  //! destructor
  virtual ~PHG4OuterHcalSubsystem( void )
  {}

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

  //! Print info (from SubsysReco)
  void Print(const std::string &what = "ALL") const;

  //! accessors (reimplemented)
  PHG4Detector* GetDetector( void ) const;
  PHG4SteppingAction* GetSteppingAction( void ) const;

  void SetLightCorrection(const double inner_radius, const double inner_corr,const double outer_radius, const double outer_corr);

  private:

  void SetDefaultParameters();

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4OuterHcalDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* steppingAction_;

};

#endif
