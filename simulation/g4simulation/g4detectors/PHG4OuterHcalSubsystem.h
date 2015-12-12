#ifndef PHG4OuterHcalSubsystem_h
#define PHG4OuterHcalSubsystem_h

#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

#include <map>
#include <string>

class PHG4OuterHcalDetector;
class PHG4OuterHcalParameters;
class PHG4Parameters;
class PHG4OuterHcalSteppingAction;
class PHG4EventAction;

class PHG4OuterHcalSubsystem: public PHG4Subsystem
{

  public:

  enum FILE_TYPE {none = 0, xml = 1, root = 2};

  //! constructor
  PHG4OuterHcalSubsystem( const std::string &name = "HCALOUT", const int layer = 0 );

  //! destructor
  virtual ~PHG4OuterHcalSubsystem( void )
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
  PHG4Detector* GetDetector( void ) const;
  PHG4SteppingAction* GetSteppingAction( void ) const;
  PHG4EventAction* GetEventAction() const {return eventAction_;}

  void SetPlaceZ(const G4double dbl);
  void SetPlace(const G4double place_x, const G4double place_y, const G4double place_z);
  void SetXRot(const G4double dbl);
  void SetYRot(const G4double dbl);
  void SetZRot(const G4double dbl);
  void SetMaterial(const std::string &mat);
  void SetActive(const int i = 1);
  void SetAbsorberActive(const int i = 1);
  void SuperDetector(const std::string &name);
  const std::string SuperDetector() {return superdetector;}

  void BlackHole(const int i=1);

  void SetTiltViaNcross(const int ncross);
  void SetTiltAngle(const double tilt);
  double GetTiltAngle() const;
  void SetInnerRadius(const double inner);
  double GetInnerRadius() const;
  void SetOuterRadius(const double outer);
  double GetOuterRadius() const;
  void SetLength(const double len);
  void SetGapWidth(const double gap);
  void SetNumScintiPlates(const int nplates);
  void SetNumScintiTiles(const int ntiles);
  void SetScintiThickness(const double thick);
  void SetScintiGap(const double scgap);
  void SetStepLimits(const double slim);


  void SetLightCorrection(const float inner_radius, const float inner_corr,
			  const float outer_radius, const float outer_corr);
  void SetLightScintModel(const bool b = true);

  void EnableFieldChecker(const int i=1) {enable_field_checker = i;}

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
  /*! defives from PHG4Detector */
  PHG4OuterHcalDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4OuterHcalSteppingAction* steppingAction_;

  //! begin/end of event action
  /*! derives from PHG4EventAction */
  PHG4EventAction *eventAction_;

  PHG4Parameters *params;

  PHG4OuterHcalParameters *paramsold;
  int enable_field_checker;
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
