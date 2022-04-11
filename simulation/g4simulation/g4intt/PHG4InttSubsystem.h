// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4INTT_PHG4INTTSUBSYSTEM_H
#define G4INTT_PHG4INTTSUBSYSTEM_H

#include <g4detectors/PHG4DetectorGroupSubsystem.h>

#include <string>   // for string
#include <utility>  // for pair
#include <vector>

class PHCompositeNode;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4InttDetector;
class PHG4SteppingAction;

class PHG4InttSubsystem : public PHG4DetectorGroupSubsystem
{
 public:
  typedef std::vector<std::pair<int, int>> vpair;

  //! constructor
  PHG4InttSubsystem(const std::string &name = "PHG4InttSubsystem", const vpair &layerconfig = vpair(0));

  //! destructor
  ~PHG4InttSubsystem() override;

  //! init
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

  //! accessors (reimplemented)
  PHG4Detector *GetDetector(void) const override;

  PHG4SteppingAction *GetSteppingAction(void) const override { return m_SteppingAction; }

  PHG4DisplayAction *GetDisplayAction() const override { return m_DisplayAction; }

  void Print(const std::string &what = "ALL") const override;

 private:
  void SetDefaultParameters() override;

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4InttDetector *m_Detector = nullptr;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction *m_SteppingAction = nullptr;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction *m_DisplayAction = nullptr;

  vpair m_LayerConfigVector;
  std::string m_DetectorType;

  std::string m_HitNodeName;
  std::string m_AbsorberNodeName;
};

#endif
