#ifndef PHG4InnerHcalSubsystem_h
#define PHG4InnerHcalSubsystem_h

#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

#include <map>
#include <set>
#include <string>

class PHG4InnerHcalDetector;
class PHG4Parameters;
class PHG4InnerHcalSteppingAction;
class PHG4EventAction;

class PHG4InnerHcalSubsystem: public PHG4Subsystem
{

  public:

  enum FILE_TYPE {xml = 1, root = 2};

  //! constructor
  PHG4InnerHcalSubsystem( const std::string &name = "HCALIN", const int layer = 0 );

  //! destructor
  virtual ~PHG4InnerHcalSubsystem( void )
  {}

  //! init
  int Init(PHCompositeNode *);

  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRun(PHCompositeNode *);

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode *);

  //! Print info (from SubsysReco)
  void Print(const std::string &what = "ALL") const;

  //! accessors (reimplemented)
  virtual PHG4Detector* GetDetector( void ) const;
  virtual PHG4SteppingAction* GetSteppingAction( void ) const;

  PHG4EventAction* GetEventAction() const {return eventAction_;}
  void SetActive(const int i = 1);
  void SetAbsorberActive(const int i = 1);
  void SuperDetector(const std::string &name);
  const std::string SuperDetector() {return superdetector;}

  void BlackHole(const int i=1);
  void SetLightCorrection(const double inner_radius, const double inner_corr,const double outer_radius, const double outer_corr);
  void set_double_param(const std::string &name, const double dval);
  void set_int_param(const std::string &name, const int ival);
  void set_string_param(const std::string &name, const std::string &sval);
  void SetDefaultParameters();
  void UpdateParametersWithMacro();
  void SaveParamsToDB();
  void SaveParamsToFile(const FILE_TYPE ftyp);

  protected:

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4InnerHcalDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingAction */
  PHG4InnerHcalSteppingAction* steppingAction_;

  //! particle tracking "stepping" action
  /*! derives from PHG4EventAction */
  PHG4EventAction *eventAction_;

  PHG4Parameters *params;

  int layer;
  std::string detector_type;
  std::string superdetector;
  std::map<const std::string, double> dparams;
  std::map<const std::string, int> iparams;
  std::map<const std::string, std::string> cparams;

};

#endif
