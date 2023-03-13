// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4OUTERHCALSUBSYSTEM_H
#define G4DETECTORS_PHG4OUTERHCALSUBSYSTEM_H

#include "PHG4DetectorSubsystem.h"

#include <string>

class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4OuterHcalDetector;
class PHG4SteppingAction;

class PHG4OuterHcalSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4OuterHcalSubsystem(const std::string& name = "HCALOUT", const int layer = 0);

  //! destructor
  ~PHG4OuterHcalSubsystem() override;

  /*!
  creates the Detector object. Creates the stepping action
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode*) override;

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode*) override;

  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const override;

  //! accessors (reimplemented)
  PHG4Detector* GetDetector(void) const override;
  PHG4SteppingAction* GetSteppingAction() const override { return m_SteppingAction; }
  PHG4DisplayAction* GetDisplayAction() const override { return m_DisplayAction; }

  void SetLightCorrection(const double inner_radius, const double inner_corr, const double outer_radius, const double outer_corr);

 private:
  void SetDefaultParameters() override;

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4OuterHcalDetector* m_Detector = nullptr;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_SteppingAction = nullptr;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction = nullptr;
};

#endif
