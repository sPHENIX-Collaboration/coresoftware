// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.1 $$
 * \date $$Date: 2014/03/24 01:36:44 $$
 */

#ifndef G4DETECTORS_PHG4SPACALSTEPPINGACTION_H
#define G4DETECTORS_PHG4SPACALSTEPPINGACTION_H

#include "LightCollectionModel.h"

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class LightCollectionModel;
class PHCompositeNode;
class PHG4CylinderCellGeomContainer;
class PHG4CylinderCellGeom_Spacalv1;
class PHG4CylinderGeomContainer;
class PHG4CylinderGeom_Spacalv3;
class PHG4SpacalDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class PHParameters;
class TowerInfoContainer;

class PHG4SpacalSteppingAction : public PHG4SteppingAction
{
 public:
  //! ctor
  explicit PHG4SpacalSteppingAction(PHG4SpacalDetector *, const PHParameters *parameters);

  //! dtor
  ~PHG4SpacalSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step *, bool) override;

  int InitWithNode(PHCompositeNode *topNode) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode *) override;

  double get_zmin() const;

  double get_zmax() const;

  void SetHitNodeName(const std::string &type, const std::string &name) override;

  int SetUpGeomNode(PHCompositeNode *topNode);

  LightCollectionModel &get_light_collection_model() { return light_collection_model; }

 private:
  bool NoHitSteppingAction(const G4Step *aStep);
  //! pointer to the detector
  PHG4SpacalDetector *m_Detector = nullptr;

  //! pointer to hit container
  PHG4HitContainer *m_HitContainer = nullptr;
  PHG4HitContainer *m_AbsorberHitContainer = nullptr;
  PHG4Hit *m_Hit = nullptr;
  PHG4HitContainer *m_CurrentHitContainer = nullptr;
  PHG4Shower *m_CurrentShower = nullptr;
  const PHParameters *m_Params = nullptr;
  int m_SaveTrackid = -1;
  int m_SavePostStepStatus = -1;
  bool m_doG4Hit = true;
  bool m_geomsetup = false;
  double m_tmin = -20.;
  double m_tmax = 60.;
  double m_dt = 100.;

  std::string m_AbsorberNodeName;
  std::string m_HitNodeName;
  std::string detector;
  std::string geonodename;
  std::string seggeonodename;

  TowerInfoContainer *m_CaloInfoContainer = nullptr;

  PHG4CylinderCellGeomContainer *_seggeo = nullptr;

  PHG4CylinderGeomContainer *_layergeo = nullptr;

  PHG4CylinderCellGeom_Spacalv1 *_geo = nullptr;

  const PHG4CylinderGeom_Spacalv3 *_layergeom = nullptr;

  LightCollectionModel light_collection_model;
};

#endif  // PHG4VHcalSteppingAction_h
