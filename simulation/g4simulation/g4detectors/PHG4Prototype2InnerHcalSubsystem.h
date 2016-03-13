#ifndef PHG4Prototype2InnerHcalSubsystem_h
#define PHG4Prototype2InnerHcalSubsystem_h

#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

#include <map>
#include <set>
#include <string>

class PHG4Prototype2InnerHcalDetector;
class PHG4Parameters;
class PHG4Prototype2InnerHcalSteppingAction;
class PHG4EventAction;

class PHG4Prototype2InnerHcalSubsystem: public PHG4Subsystem
{

  public:

  enum FILE_TYPE {none = 0, xml = 1, root = 2};

  //! constructor
  PHG4Prototype2InnerHcalSubsystem( const std::string &name = "HCALIN", const int layer = 0 );

  //! destructor
  virtual ~PHG4Prototype2InnerHcalSubsystem( void )
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
  void SetAbsorberTruth(const int i = 1);
  void SuperDetector(const std::string &name);
  const std::string SuperDetector() {return superdetector;}

  void BlackHole(const int i=1);
  void SetLightCorrection(const double inner_radius, const double inner_corr,const double outer_radius, const double outer_corr);
  void set_double_param(const std::string &name, const double dval);
  double get_double_param(const std::string &name) const;
  void set_int_param(const std::string &name, const int ival);
  int get_int_param(const std::string &name) const;
  void set_string_param(const std::string &name, const std::string &sval);
  std::string get_string_param(const std::string &name) const;
  void SetDefaultParameters();
  void UpdateParametersWithMacro();
  void UseDB(const int i = 1) {usedb = i;}
  void UseCalibFiles(const FILE_TYPE ftyp) {filetype = ftyp;}
  int SaveParamsToDB();
  int ReadParamsFromDB();
  int SaveParamsToFile(const FILE_TYPE ftyp);
  int ReadParamsFromFile(const FILE_TYPE ftyp);
  void SetCalibrationFileDir(const std::string &calibdir) {calibfiledir = calibdir;}

  protected:

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4Prototype2InnerHcalDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingAction */
  PHG4Prototype2InnerHcalSteppingAction* steppingAction_;

  //! particle tracking "stepping" action
  /*! derives from PHG4EventAction */
  PHG4EventAction *eventAction_;

  PHG4Parameters *params;

  int layer;

  int usedb;
  FILE_TYPE filetype;
  std::string detector_type;
  std::string superdetector;
  std::string calibfiledir;
  std::map<const std::string, double> dparams;
  std::map<const std::string, int> iparams;
  std::map<const std::string, std::string> cparams;
  std::map<const std::string, double> default_double;
  std::map<const std::string, int> default_int;
  std::map<const std::string, std::string> default_string;

};

#endif
