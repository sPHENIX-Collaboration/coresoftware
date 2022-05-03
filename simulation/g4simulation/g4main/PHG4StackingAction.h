// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4STACKINGACTION_H
#define G4MAIN_PHG4STACKINGACTION_H

#include <Geant4/G4UserStackingAction.hh>

#include <set>
#include <string>


class G4Track;
class PHCompositeNode;
class PHG4Hit;

class PHG4StackingAction
{
 public:
  PHG4StackingAction(const std::string& name, const int i = 0);
  virtual ~PHG4StackingAction()
  {
  }

  //! stacking action. This gets called by the stacking action when a new track is generated
  // It can classify those tracks (urgent, wait, kill, postpone to next event)
  // default return in G4 is fUrgent
  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* /*aTrack*/) {return fUrgent;}
  virtual void PrepareNewEvent() {}
  virtual void Verbosity(const int i) { m_Verbosity = i; }
  virtual int Verbosity() const { return m_Verbosity; }
  virtual int Init() { return 0; }

  virtual void SetInterfacePointers(PHCompositeNode*) { return; }
  virtual void Print(const std::string& /*what*/) const { return; }
  std::string GetName() const { return m_Name; }
  void SetName(const std::string& name) { m_Name = name; }

 private:
  int m_Verbosity = 0;
  std::string m_Name;
};

#endif  // G4MAIN_PHG4STACKINGACTION_H
