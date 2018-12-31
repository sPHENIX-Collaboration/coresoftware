#ifndef G4MAIN_PHG4RECO_H
#define G4MAIN_PHG4RECO_H

#include <g4decayer/EDecayType.hh>

#include <fun4all/SubsysReco.h>

#include <phfield/PHFieldConfig.h>

#include <list>

// Forward declerations
class PHCompositeNode;
class G4RunManager;
class PHG4PrimaryGeneratorAction;
class PHG4PhenixDetector;
class PHG4PhenixEventAction;
class PHG4PhenixSteppingAction;
class PHG4PhenixTrackingAction;
class PHG4Subsystem;
class PHG4EventGenerator;
class G4VModularPhysicsList;
class G4TBMagneticFieldSetup;
class G4VUserPrimaryGeneratorAction;
class PHG4UIsession;

// for the G4 cmd interface and the graphics
class G4UImanager;
class G4VisManager;

/*!
  \class   PHG4Reco
  \ingroup supermodules
  \brief   Runs G4 as a subsystem
*/
class PHG4Reco : public SubsysReco
{
 public:
  //! constructor
  PHG4Reco(const std::string &name = "PHG4RECO");

  //! destructor
  virtual ~PHG4Reco();

  //! full initialization
  int Init(PHCompositeNode *);

  int InitRun(PHCompositeNode *topNode);

  //! event processing method
  int process_event(PHCompositeNode *);

  //! Clean up after each event.
  int ResetEvent(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

  //! print info
  void Print(const std::string &what = std::string()) const;

  //! register subsystem
  void registerSubsystem(PHG4Subsystem *subsystem)
  {
    subsystems_.push_back(subsystem);
  }

  //! interface to G4 cmd interpreter
  int ApplyCommand(const std::string &cmd);

  //! start the gui
  int StartGui();

  int InitField(PHCompositeNode *topNode);

  //! set default magnetic field strength with a constant magnetic field. Only valid if set_field_map() is not used. If available, Field map setting on DST take higher priority.
  void set_field(const float tesla)
  {
    magfield = tesla;
  }

  //! Set default field map. If available, Field map setting on DST take higher priority.
  //! \param[in] fmap  Field map ROOT file
  //! \param[in] dim   Field map format. See PHFieldConfig::FieldConfigTypes for available formats.
  void set_field_map(const std::string &fmap, const PHFieldConfig::FieldConfigTypes dim)
  {
    fieldmapfile = fmap;
    mapdim = dim;
  }

  //! set default scaling factor for input magnetic field map. If available, Field map setting on DST take higher priority.
  void set_field_rescale(const float rescale) { magfield_rescale = rescale; }

  void set_decayer_active(bool b) { active_decayer_ = b; }
  void set_force_decay(EDecayType force_decay_type)
  {
    active_decayer_ = true;
    active_force_decay_ = true;
    force_decay_type_ = force_decay_type;
  }

  //! Save geometry from Geant4 to DST
  void save_DST_geometry(bool b) { save_DST_geometry_ = b; }
  void SetWorldSizeX(const double sx) { WorldSize[0] = sx; }
  void SetWorldSizeY(const double sy) { WorldSize[1] = sy; }
  void SetWorldSizeZ(const double sz) { WorldSize[2] = sz; }
  double GetWorldSizeX() const { return WorldSize[0]; }
  double GetWorldSizeY() const { return WorldSize[1]; }
  double GetWorldSizeZ() const { return WorldSize[2]; }
  void SetWorldShape(const std::string &s) { worldshape = s; }
  void SetWorldMaterial(const std::string &s) { worldmaterial = s; }
  void SetPhysicsList(const std::string &s) { physicslist = s; }
  void set_rapidity_coverage(const double eta);

  int setupInputEventNodeReader(PHCompositeNode *);

  static void G4Seed(const unsigned int i);

  // this is an ugly hack to get Au ions working for CAD
  // our particle generators have pdg build in which doesn't work
  // with ions, so the generator action has to be replaced
  // which is hardcoded in PHG4Reco (it has to be created after
  // the physics lists are instantiated
  void setGeneratorAction(G4VUserPrimaryGeneratorAction *action);

  PHG4Subsystem *getSubsystem(const std::string &name);

  void Dump_GDML(const std::string &filename);

  void G4Verbosity(const int i);

  //! disable event/track/stepping actions to reduce resource consumption for G4 running only. E.g. dose analysis
  void setDisableUserActions(bool b = true) {m_disableUserActions = b;}

 protected:
  int InitUImanager();
  void DefineMaterials();
  void DefineRegions();
  float magfield;
  float magfield_rescale;
  double WorldSize[3];

  //! magnetic field
  G4TBMagneticFieldSetup *field_;

  //! pointer to geant run manager
  G4RunManager *runManager_;

  //! pointer to geant ui session
  PHG4UIsession *uisession_;

  //! pointer to detector
  PHG4PhenixDetector *detector_;

  //! pointer to main event action
  PHG4PhenixEventAction *eventAction_;

  //! pointer to main stepping action
  PHG4PhenixSteppingAction *steppingAction_;

  //! pointer to main tracking action
  PHG4PhenixTrackingAction *trackingAction_;

  //! event generator (read from PHG4INEVENT node)
  PHG4PrimaryGeneratorAction *generatorAction_;

  //! list of subsystems
  typedef std::list<PHG4Subsystem *> SubsystemList;
  SubsystemList subsystems_;

  // visualization
  G4VisManager *visManager;

  double _eta_coverage;
  PHFieldConfig::FieldConfigTypes mapdim;
  std::string fieldmapfile;
  std::string worldshape;
  std::string worldmaterial;
  std::string physicslist;

  // settings for the external Pythia6 decayer
  bool active_decayer_;          //< turn on/off decayer
  bool active_force_decay_;      //< turn on/off force decay channels
  EDecayType force_decay_type_;  //< forced decay channel setting

  bool save_DST_geometry_;
  bool m_disableUserActions;

};

#endif
