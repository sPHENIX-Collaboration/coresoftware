// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4INNERHCALSUBSYSTEM_H
#define G4DETECTORS_PHG4INNERHCALSUBSYSTEM_H

#include "PHG4DetectorSubsystem.h"

#include <string>

class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4InnerHcalDetector;
class PHG4SteppingAction;

class PHG4InnerHcalSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4InnerHcalSubsystem(const std::string& name = "HCALIN", const int layer = 0);

  //! destructor
  ~PHG4InnerHcalSubsystem() override;

  /*!
  creates the m_Detector object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
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
  PHG4SteppingAction* GetSteppingAction(void) const override { return m_SteppingAction; }
  PHG4DisplayAction* GetDisplayAction() const override { return m_DisplayAction; }

  void SetLightCorrection(const double inner_radius, const double inner_corr, const double outer_radius, const double outer_corr);

 private:
  void SetDefaultParameters() override;

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4InnerHcalDetector* m_Detector = nullptr;

  //! detector "stepping" action, executes after every G4 step
  /*! derives from PHG4SteppingAction */
  PHG4SteppingAction* m_SteppingAction = nullptr;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction = nullptr;
};

#endif  // G4DETECTORS_PHG4INNERHCALSUBSYSTEM_H
