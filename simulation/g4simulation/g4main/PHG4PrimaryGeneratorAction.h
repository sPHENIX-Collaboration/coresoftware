// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PRIMARYGENERATORACTION_H
#define G4MAIN_PHG4PRIMARYGENERATORACTION_H

#include <Geant4/G4VUserPrimaryGeneratorAction.hh>

class PHG4InEvent;

class PHG4PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
  PHG4PrimaryGeneratorAction()
    : verbosity(0)
    , inEvent(0)
  {
  }

  virtual ~PHG4PrimaryGeneratorAction()
  {
  }

  virtual void GeneratePrimaries(G4Event* anEvent);

  //! set top node (from where particle list is retrieved for passing to geant
  void SetInEvent(PHG4InEvent* const inevt)
  {
    inEvent = inevt;
  }

  //! Set/Get verbosity
  void Verbosity(const int val) { verbosity = val; }
  int Verbosity() const { return verbosity; }

 protected:
  int verbosity;

 private:
  //! temporary pointer to input event on node tree
  PHG4InEvent* inEvent;
};

#endif  // PHG4PrimaryGeneratorAction_H__
