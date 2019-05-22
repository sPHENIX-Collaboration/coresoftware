// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4INTT_PHG4INTTSUBSYSTEM_H
#define G4INTT_PHG4INTTSUBSYSTEM_H

#include <g4detectors/PHG4DetectorGroupSubsystem.h>

#include <vector>

class PHG4INTTDetector;
class PHG4INTTSteppingAction;
class PHG4DisplayAction;

class PHG4INTTSubsystem : public PHG4DetectorGroupSubsystem
{
 public:
  typedef std::vector<std::pair<int, int>> vpair;

  //! constructor
  PHG4INTTSubsystem(const std::string &name = "INTT", const vpair &layerconfig = vpair(0));

  //! destructor
  virtual ~PHG4INTTSubsystem();

  //! init
  /*!
  called during InitRun (the original InitRun does common setup and calls this one)
  creates the detector object 
  ceates the stepping action 
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

  PHG4SteppingAction *GetSteppingAction(void) const { return m_SteppingAction; }

  PHG4DisplayAction* GetDisplayAction() const { return m_DisplayAction; }

  void Print(const std::string &what = "ALL") const;

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4INTTDetector *m_Detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction *m_SteppingAction;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction* m_DisplayAction;

  vpair m_LayerConfigVector;
  std::string m_DetectorType;
};

#endif
