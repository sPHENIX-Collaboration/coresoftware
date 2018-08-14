#ifndef PHG4SiliconTrackerSubsystem_h
#define PHG4SiliconTrackerSubsystem_h

#include "PHG4DetectorGroupSubsystem.h"

#include <Geant4/G4String.hh>
#include <Geant4/G4Types.hh>

#include <utility>
#include <vector>

class PHG4SiliconTrackerDetector;
class PHG4SiliconTrackerSteppingAction;

class PHG4SiliconTrackerSubsystem : public PHG4DetectorGroupSubsystem
{
 public:
  typedef std::vector<std::pair<int, int>> vpair;

  //! constructor
  PHG4SiliconTrackerSubsystem(const std::string &name = "SILICONTRACKER", const vpair &layerconfig = vpair(0));  

  PHG4SiliconTrackerSubsystem(const double sensor_radius_inner_[], const double sensor_radius_outer_[], const std::string &name = "SILICONTRACKER", const vpair &layerconfig = vpair(0));

  //! destructor
  virtual ~PHG4SiliconTrackerSubsystem(void)
  {
  }

  //! init
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
  PHG4Detector *GetDetector(void) const;
  PHG4SteppingAction *GetSteppingAction(void) const { return steppingAction_; }
  void Print(const std::string &what = "ALL") const;


 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4SiliconTrackerDetector *detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction *steppingAction_;

  std::vector<std::pair<int, int>> layerconfig_;
  std::string detector_type;

  int nlayers;

  double sensor_radius_inner[4];
  double sensor_radius_outer[4];

};

#endif
