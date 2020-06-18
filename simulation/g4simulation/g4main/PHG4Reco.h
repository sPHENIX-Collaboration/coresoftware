// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4RECO_H
#define G4MAIN_PHG4RECO_H

#include <g4decayer/EDecayType.hh>

#include <fun4all/SubsysReco.h>

#include <phfield/PHFieldConfig.h>

#include <list>
#include <string>  // for string

// Forward declerations
class G4RunManager;
class G4TBMagneticFieldSetup;
class G4UImanager;
class G4UImessenger;
class G4VisManager;
class PHCompositeNode;
class PHG4DisplayAction;
class PHG4PhenixDetector;
class PHG4PhenixEventAction;
class PHG4PhenixSteppingAction;
class PHG4PhenixTrackingAction;
class PHG4PrimaryGeneratorAction;
class PHG4Subsystem;
class PHG4UIsession;

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

  //! print info
  void Print(const std::string &what = std::string()) const;

  //! register subsystem
  void registerSubsystem(PHG4Subsystem *subsystem)
  {
    m_SubsystemList.push_back(subsystem);
  }

  //! interface to G4 cmd interpreter
  int ApplyCommand(const std::string &cmd);

  //! start the gui
  int StartGui();

  int InitField(PHCompositeNode *topNode);

  //! set default magnetic field strength with a constant magnetic field. Only valid if set_field_map() is not used. If available, Field map setting on DST take higher priority.
  void set_field(const float tesla)
  {
    m_MagneticField = tesla;
  }

  //! Set default field map. If available, Field map setting on DST take higher priority.
  //! \param[in] fmap  Field map ROOT file
  //! \param[in] dim   Field map format. See PHFieldConfig::FieldConfigTypes for available formats.
  void set_field_map(const std::string &fmap, const PHFieldConfig::FieldConfigTypes dim)
  {
    m_FieldMapFile = fmap;
    m_FieldConfigType = dim;
  }

  //! set default scaling factor for input magnetic field map. If available, Field map setting on DST take higher priority.
  void set_field_rescale(const float rescale) { m_MagneticFieldRescale = rescale; }

  void set_decayer_active(bool b) { m_ActiveDecayerFlag = b; }
  void set_force_decay(EDecayType force_decay_type)
  {
    m_ActiveDecayerFlag = true;
    m_ActiveForceDecayFlag = true;
    m_ForceDecayType = force_decay_type;
  }

  //! Save geometry from Geant4 to DST
  void save_DST_geometry(bool b) { m_SaveDstGeometryFlag = b; }
  void SetWorldSizeX(const double sx) { m_WorldSize[0] = sx; }
  void SetWorldSizeY(const double sy) { m_WorldSize[1] = sy; }
  void SetWorldSizeZ(const double sz) { m_WorldSize[2] = sz; }
  double GetWorldSizeX() const { return m_WorldSize[0]; }
  double GetWorldSizeY() const { return m_WorldSize[1]; }
  double GetWorldSizeZ() const { return m_WorldSize[2]; }
  void SetWorldShape(const std::string &s) { m_WorldShape = s; }
  void SetWorldMaterial(const std::string &s) { m_WorldMaterial = s; }
  void SetPhysicsList(const std::string &s) { m_PhysicsList = s; }
  void set_rapidity_coverage(const double eta);

  int setupInputEventNodeReader(PHCompositeNode *);

  static void G4Seed(const unsigned int i);

  PHG4Subsystem *getSubsystem(const std::string &name);
  PHG4DisplayAction *GetDisplayAction() { return m_DisplayAction; }
  void Dump_GDML(const std::string &filename);

  void G4Verbosity(const int i);

  //! disable event/track/stepping actions to reduce resource consumption for G4 running only. E.g. dose analysis
  void setDisableUserActions(bool b = true) { m_disableUserActions = b; }
  void ApplyDisplayAction();

 private:
  static void g4guithread(void *ptr);
  int InitUImanager();
  void DefineMaterials();
  void DefineRegions();

  float m_MagneticField;
  float m_MagneticFieldRescale;
  double m_WorldSize[3];

  //! magnetic field
  G4TBMagneticFieldSetup *m_Field;

  //! pointer to geant run manager
  G4RunManager *m_RunManager;

  //! pointer to geant ui session
  PHG4UIsession *m_UISession;

  //! pointer to detector
  PHG4PhenixDetector *m_Detector;

  //! pointer to main event action
  PHG4PhenixEventAction *m_EventAction;

  //! pointer to main stepping action
  PHG4PhenixSteppingAction *m_SteppingAction;

  //! pointer to main tracking action
  PHG4PhenixTrackingAction *m_TrackingAction;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction *m_DisplayAction;

  //! event generator (read from PHG4INEVENT node)
  PHG4PrimaryGeneratorAction *m_GeneratorAction;

  //! list of subsystems
  std::list<PHG4Subsystem *> m_SubsystemList;

  // visualization
  G4VisManager *m_VisManager;

  // Message interface to Fun4All
  G4UImessenger *m_Fun4AllMessenger;

  // for the G4 cmd line interface
  G4UImanager *m_UImanager;
  double m_EtaCoverage;
  PHFieldConfig::FieldConfigTypes m_FieldConfigType;
  std::string m_FieldMapFile;
  std::string m_WorldShape;
  std::string m_WorldMaterial;
  std::string m_PhysicsList;

  // settings for the external Pythia6 decayer
  bool m_ActiveDecayerFlag;     //< turn on/off decayer
  bool m_ActiveForceDecayFlag;  //< turn on/off force decay channels
  EDecayType m_ForceDecayType;  //< forced decay channel setting

  bool m_SaveDstGeometryFlag;
  bool m_disableUserActions;
};

#endif
