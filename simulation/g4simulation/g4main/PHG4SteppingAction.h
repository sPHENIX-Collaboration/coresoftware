#ifndef PHG4SteppingAction_h
#define PHG4SteppingAction_h

#include <map>
#include <set>
#include <string>

class G4Step;
class PHCompositeNode;
class PHG4Hit;

class PHG4SteppingAction
{

  public:
  PHG4SteppingAction( const int i = 0 ):
    verbosity(i),
    name("NONAME")
  {}

  virtual ~PHG4SteppingAction()
  {}

  //! stepping action. This defines processing of a single step in a given volume
  /*!
  returns true if hits associated to this step was processed by this detector
  \param step pointer to the geant 4 step class
  \param was_used: true if the hit was already used by a previously registered subsystem
  */
  virtual bool UserSteppingAction(const G4Step* step, bool was_used ) = 0;

  virtual void Verbosity(const int i) {verbosity = i;}
  int Verbosity() const {return verbosity;}

  virtual int Init() {return 0;}

  //! get scintillation photon count. It require a custom set SCINTILLATIONYIELD property to work
  virtual double GetScintLightYield(const G4Step* step);

  //! get amount of energy that can make scintillation light, in Unit of GeV.
  virtual double GetVisibleEnergyDeposition(const G4Step* step);

  //! Extract local coordinate of the hit and save to PHG4Hit
  virtual void StoreLocalCoordinate(PHG4Hit * hit, const G4Step* step, const bool do_prepoint, const bool do_postpoint);

  virtual void SetInterfacePointers( PHCompositeNode* ) {return;}

  virtual void Print(const std::string &what) const {return;}

  void SetOpt(const std::string &name, const int i) {opt_int[name] = i;}
  bool IntOptExist(const std::string &name);
  int GetIntOpt(const std::string &name);
  std::string GetName() const {return name;}
 protected:
  int verbosity;
  std::string name;

 private:
  std::set<std::string> _ScintLightYieldMissingMaterial;
  std::map<std::string, int> opt_int;

};


#endif
