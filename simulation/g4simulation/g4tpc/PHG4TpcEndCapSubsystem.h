// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHG4TPCENDCAPSUBSYSTEM_H
#define PHG4TPCENDCAPSUBSYSTEM_H

#include <g4detectors/PHG4DetectorSubsystem.h>

#include <string>  // for allocator, string

class PHCompositeNode;
class PHG4Detector;
class PHG4TpcEndCapDetector;
class PHG4SteppingAction;
class PHG4DisplayAction;

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
  ~PHG4TpcEndCapSubsystem() override;

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

  PHG4DisplayAction* GetDisplayAction() const override { return m_DisplayAction; }

 protected:
  // \brief Set default parameter values
  void SetDefaultParameters() override;

 private:
  //! detector construction
  /*! derives from PHG4Detector */
  PHG4TpcEndCapDetector* m_Detector = nullptr;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_SteppingAction = nullptr;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction = nullptr;

  //! node name for the PHG4Hits
  std::string m_HitNodeName;
};

#endif  // PHG4TPCENDCAPSUBSYSTEM_H
