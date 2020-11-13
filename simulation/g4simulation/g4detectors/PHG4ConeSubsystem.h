// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CONESUBSYSTEM_H
#define G4DETECTORS_PHG4CONESUBSYSTEM_H

#include "PHG4DetectorSubsystem.h"

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

#include <string>                  // for string

class PHCompositeNode;
class PHG4ConeDetector;
class PHG4ConeSteppingAction;
class PHG4Detector;
class PHG4SteppingAction;

class PHG4ConeSubsystem: public PHG4DetectorSubsystem
{

  public:

  //! constructor
  PHG4ConeSubsystem( const std::string &name = "CONE", const int layer = 0 );

  //! destructor
  virtual ~PHG4ConeSubsystem( void )
  {}

  //! init runwise stuff
  /*!
  creates the m_Detector object
  creates the stepping action
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
  PHG4Detector* GetDetector( void ) const override;
  PHG4SteppingAction* GetSteppingAction( void ) const override { return m_SteppingAction; };

  //!set inner and outter radius1
  void SetR1(const G4double min, const G4double max);

  //!set inner and outter radius2
  void SetR2(const G4double min, const G4double max);

  //! set length in Z
  void SetZlength(const G4double a);

  //! set phi offset and extention
  void SetPhi(const G4double a, const G4double b);

  //! set rmaximum and minimums according to the eta range 
  void Set_eta_range(G4double etaMin, G4double etaMax);

  void SetPlaceZ(const G4double dbl);
  void SetPlace(const G4double place_x, const G4double place_y, const G4double place_z);

  void SetZRot(const G4double dbl);
  void SetMaterial(const std::string &mat);

// this method is used to check if it can be used as mothervolume
// Subsystems which can be mothervolume need to implement this 
// and return true
  virtual bool CanBeMotherSubsystem() const  override {return true;}

  private:

  void SetDefaultParameters() override;

  //! detector geometry
  /*! derives from PHG4Detector */
  PHG4ConeDetector* m_Detector = nullptr;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* m_SteppingAction = nullptr;

//  int layer;
};

#endif
