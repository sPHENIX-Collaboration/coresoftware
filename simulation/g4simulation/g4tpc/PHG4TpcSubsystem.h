// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4TPC_PHG4TPCSUBSYSTEM_H
#define G4TPC_PHG4TPCSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

#include <string>

class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4SteppingAction;
class PHG4TpcDetector;

class PHG4TpcSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4TpcSubsystem(const std::string &name = "TPC", const int layer = 0);

  //! destructor
  ~PHG4TpcSubsystem() override;

  /*!
  called during InitRun (the original InitRun does common setup and calls this one)
  creates the detector object 
  ceates the stepping action 
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode *) override;

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode *) override;

  //! Print info (from SubsysReco)
  void Print(const std::string &what = "ALL") const override;

  //! accessors (reimplemented)
  PHG4Detector *GetDetector(void) const override;

  PHG4SteppingAction *GetSteppingAction(void) const override { return m_SteppingAction; }

  PHG4DisplayAction *GetDisplayAction() const override { return m_DisplayAction; }

 private:
  void SetDefaultParameters() override;

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4TpcDetector *m_Detector = nullptr;

  //! detector "stepping" action, executes after every G4 step
  /*! derives from PHG4SteppingAction */
  PHG4SteppingAction *m_SteppingAction = nullptr;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction *m_DisplayAction = nullptr;

  std::string m_HitNodeName;
  std::string m_AbsorberNodeName;
};

#endif
