// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MVTX_PHG4MVTXSUBSYSTEM_H
#define G4MVTX_PHG4MVTXSUBSYSTEM_H

#include <g4detectors/PHG4DetectorGroupSubsystem.h>

#include <cmath>   // for asin
#include <string>  // for string

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
  ~PHG4MvtxSubsystem() override;

  //! InitRunSubsystem
  /*!
  called during InitRun (the original InitRun does common setup and calls this one)
  creates the detector object
  creates the stepping action
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode*) override;

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode*) override;

  //! accessors (reimplemented)
  PHG4Detector* GetDetector(void) const override;
  PHG4SteppingAction* GetSteppingAction(void) const override;

  PHG4DisplayAction* GetDisplayAction() const override { return m_DisplayAction; }

 private:
  void SetDefaultParameters() override;
  static double radii2Turbo(double rMin, double rMid, double rMax, double sensW)
  {
    // compute turbo angle from radii and sensor width
    return std::asin((rMax * rMax - rMin * rMin) / (2 * rMid * sensW));
  }
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
