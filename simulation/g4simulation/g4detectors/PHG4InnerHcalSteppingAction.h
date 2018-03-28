#ifndef PHG4VInnerHcalSteppingAction_h
#define PHG4VInnerHcalSteppingAction_h

#include <g4main/PHG4SteppingAction.h>

class G4VPhysicalVolume;
class PHG4InnerHcalDetector;
class PHParameters;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

class PHG4InnerHcalSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4InnerHcalSteppingAction(PHG4InnerHcalDetector *, const PHParameters *parameters);

  //! destructor
  virtual ~PHG4InnerHcalSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step *, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode *);

  double GetLightCorrection(const double r) const;

 private:
  //! pointer to the detector
  PHG4InnerHcalDetector *detector_;

  //! pointer to hit container
  PHG4HitContainer *hits_;
  PHG4HitContainer *absorberhits_;
  PHG4Hit *hit;
  const PHParameters *params;
  PHG4HitContainer *savehitcontainer;
  PHG4Shower *saveshower;
  G4VPhysicalVolume *savevolpre;
  G4VPhysicalVolume *savevolpost;
  int savetrackid;
  int saveprestepstatus;
  int savepoststepstatus;
  // since getting parameters is a map search we do not want to
  // do this in every step, the parameters used are cached
  // in the following variables
  int absorbertruth;
  int IsActive;
  int IsBlackHole;
  int n_scinti_plates;
  int light_scint_model;

  double light_balance_inner_corr;
  double light_balance_inner_radius;
  double light_balance_outer_corr;
  double light_balance_outer_radius;
};

#endif  // PHG4InnerHcalSteppingAction_h
