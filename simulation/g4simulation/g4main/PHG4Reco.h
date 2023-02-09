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
class PHG4PhenixStackingAction;
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
  ~PHG4Reco() override;

  //! full initialization
  int Init(PHCompositeNode *) override;

  int InitRun(PHCompositeNode *topNode) override;

  //! event processing method
  int process_event(PHCompositeNode *) override;

  //! Clean up after each event.
  int ResetEvent(PHCompositeNode *) override;

  //! print info
  void Print(const std::string &what = std::string()) const override;

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

//  void set_decayer_active(bool b) { m_ActiveDecayerFlag = b; }
  void set_force_decay(EDecayType force_decay_type)
  {
//    m_ActiveDecayerFlag = true;
    m_ActiveForceDecayFlag = true;
    m_ForceDecayType = force_decay_type;
  }

  //! export geometry to root file
  void export_geometry( bool b, const std::string& filename = "sPHENIXGeom.root" )
  {
    m_ExportGeometry = b;
    m_ExportGeomFilename = filename;
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
  void Dump_G4_GDML(const std::string &filename);

  void G4Verbosity(const int i);

  //! disable event/track/stepping actions to reduce resource consumption for G4 running only. E.g. dose analysis
  void setDisableUserActions(bool b = true) { m_disableUserActions = b; }
  void ApplyDisplayAction();

 private:
  static void g4guithread(void *ptr);
  int InitUImanager();
  void DefineMaterials();
  void DefineRegions();

  float m_MagneticField = 0.;
  float m_MagneticFieldRescale = 1.0;
  double m_WorldSize[3];

  //! magnetic field
  G4TBMagneticFieldSetup *m_Field = nullptr;

  //! pointer to geant run manager
  G4RunManager *m_RunManager = nullptr;

  //! pointer to geant ui session
  PHG4UIsession *m_UISession = nullptr;

  //! pointer to detector
  PHG4PhenixDetector *m_Detector = nullptr;

  //! pointer to main event action
  PHG4PhenixEventAction *m_EventAction = nullptr;

  //! pointer to main stacking action
  PHG4PhenixStackingAction *m_StackingAction = nullptr;

  //! pointer to main stepping action
  PHG4PhenixSteppingAction *m_SteppingAction = nullptr;

  //! pointer to main tracking action
  PHG4PhenixTrackingAction *m_TrackingAction = nullptr;

  //! display attribute setting
  /*! derives from PHG4DisplayAction */
  PHG4DisplayAction *m_DisplayAction = nullptr;

  //! event generator (read from PHG4INEVENT node)
  PHG4PrimaryGeneratorAction *m_GeneratorAction = nullptr;

  //! list of subsystems
  std::list<PHG4Subsystem *> m_SubsystemList;

  // visualization
  G4VisManager *m_VisManager = nullptr;

  // Message interface to Fun4All
  G4UImessenger *m_Fun4AllMessenger = nullptr;

  // for the G4 cmd line interface
  G4UImanager *m_UImanager = nullptr;
  double m_EtaCoverage = 1.0;
  PHFieldConfig::FieldConfigTypes m_FieldConfigType = PHFieldConfig::kFieldUniform;
  std::string m_FieldMapFile = "NONE";
  std::string m_WorldShape = "G4Tubs";
  std::string m_WorldMaterial = "G4_AIR";
  std::string m_PhysicsList = "FTFP_BERT";

  bool m_ExportGeometry = false;
  std::string m_ExportGeomFilename = "sPHENIXGeom.root";
 
  // settings for the external Pythia6 decayer
  //bool m_ActiveDecayerFlag = true;     //< turn on/off decayer
  bool m_ActiveForceDecayFlag = false;  //< turn on/off force decay channels

  enum DecayerOptions
  {
    kGEANTInternalDecayer = 0,
    kPYTHIA6Decayer = 1,
    kEvtGenDecayer = 2,

  };  // Decayer Option for User to Choose: 0 - GEANT 4 Internal Decayer (with momentum conservation issues), 1, PYTHIA 6 Decayer, 2 - EvtGen Decayer

  DecayerOptions m_Decayer = kEvtGenDecayer;  // Here we use EvtGen as default
  EDecayType m_ForceDecayType = kAll;  //< forced decay channel setting

  bool m_SaveDstGeometryFlag = true;
  bool m_disableUserActions = false;
};

#endif
