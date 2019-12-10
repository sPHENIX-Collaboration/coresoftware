// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MVTX_PHG4MVTXSUBSYSTEM_H
#define G4MVTX_PHG4MVTXSUBSYSTEM_H

#include <g4detectors/PHG4DetectorGroupSubsystem.h>

#include <string>                                    // for string

class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4MvtxDetector;
class PHG4SteppingAction;

class PHG4MvtxSubsystem : public PHG4DetectorGroupSubsystem
{
 public:
  //! constructor
  PHG4MvtxSubsystem(const std::string& name = "PHG4MvtxSubsystem", const int _n_layers = 3);

  //! destructor
  virtual ~PHG4MvtxSubsystem();

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

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4MvtxDetector* m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* steppingAction_;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction;

  // These are passed on to the detector class
  unsigned int n_layers;

  std::string detector_type;
};

#endif
