// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MVTX_PHG4MVTXSUBSYSTEM_H
#define G4MVTX_PHG4MVTXSUBSYSTEM_H

#include <g4detectors/PHG4DetectorGroupSubsystem.h>

#include <Geant4/G4Types.hh>

class PHG4MVTXDetector;
class PHG4MVTXSteppingAction;
class PHG4DisplayAction;
class PHG4EventAction;

class PHG4MVTXSubsystem : public PHG4DetectorGroupSubsystem
{
 public:
  //! constructor
  PHG4MVTXSubsystem(const std::string& name = "MVTX", const int _n_layers = 3);

  //! destructor
  virtual ~PHG4MVTXSubsystem();

  //! InitRunSubsystem
  /*!
  called during InitRun (the original InitRun does common setup and calls this one)
  creates the detector object 
  ceates the stepping action 
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode*);

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode*);

  //! accessors (reimplemented)
  virtual PHG4Detector* GetDetector(void) const;
  virtual PHG4SteppingAction* GetSteppingAction(void) const;

  PHG4DisplayAction* GetDisplayAction() const { return m_DisplayAction; }

  PHG4EventAction* GetEventAction() const { return eventAction_; }

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4MVTXDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4MVTXSteppingAction* steppingAction_;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction;

  PHG4EventAction* eventAction_;

  // These are passed on to the detector class
  G4double layer_nominal_radius;
  unsigned int n_layers;

  std::string detector_type;
};

#endif
