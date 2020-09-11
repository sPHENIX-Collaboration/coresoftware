// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHG4TPCENDCAPSUBSYSTEM_H
#define PHG4TPCENDCAPSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

class PHCompositeNode;
class PHG4Detector;
class PHG4TpcEndCapDetector;
class PHG4SteppingAction;

/**
   * \brief Detector Subsystem module
   *
   * The detector is constructed and registered via PHG4TpcEndCapDetector
   *
   *
   * \see PHG4TpcEndCapDetector
   * \see PHG4TpcEndCapSubsystem
   *
   */
class PHG4TpcEndCapSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4TpcEndCapSubsystem(const std::string& name = "PHG4TpcEndCap");

  //! destructor
  virtual ~PHG4TpcEndCapSubsystem() {}

  /*!
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
  PHG4Detector* GetDetector() const override;

  PHG4SteppingAction* GetSteppingAction() const override { return m_SteppingAction; }
  //! Print info (from SubsysReco)
  void Print(const std::string& what = "ALL") const override;

 protected:
  // \brief Set default parameter values
  void SetDefaultParameters() override;

 private:
  //! detector construction
  /*! derives from PHG4Detector */
  PHG4TpcEndCapDetector  *m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction *m_SteppingAction;
};

#endif // PHG4TPCENDCAPSUBSYSTEM_H
