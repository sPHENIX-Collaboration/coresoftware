#ifndef G4TPC_PHG4TPCPADPLANE_H
#define G4TPC_PHG4TPCPADPLANE_H

#include <fun4all/SubsysReco.h>

#include <g4main/PHG4HitContainer.h>
#include "TpcClusterBuilder.h"

#include <phparameter/PHParameterInterface.h>

#include <string>  // for string

class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class PHCompositeNode;
class PHG4TpcCylinderGeomContainer;
class TNtuple;

class PHG4TpcPadPlane : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4TpcPadPlane(const std::string &name = "PHG4TpcPadPlane");

  ~PHG4TpcPadPlane() override {}

  int process_event(PHCompositeNode *) final
  {
    return 0;
  }
  int InitRun(PHCompositeNode *topNode) override;
  virtual void SetDriftVelocity(double /*vd*/) {return;}
  virtual int CreateReadoutGeometry(PHCompositeNode * /*topNode*/, PHG4TpcCylinderGeomContainer * /*seggeo*/) { return 0; }
  virtual void UpdateInternalParameters() { return; }
  //  virtual void MapToPadPlane(PHG4CellContainer * /*g4cells*/, const double /*x_gem*/, const double /*y_gem*/, const double /*t_gem*/, const unsigned int /*side*/, PHG4HitContainer::ConstIterator /*hiter*/, TNtuple * /*ntpad*/, TNtuple * /*nthit*/) {}
  virtual void MapToPadPlane(TpcClusterBuilder* /*builder*/, TrkrHitSetContainer * /*single_hitsetcontainer*/, TrkrHitSetContainer * /*hitsetcontainer*/, TrkrHitTruthAssoc * /*hittruthassoc*/, const double /*x_gem*/, const double /*y_gem*/, const double /*t_gem*/, const unsigned int /*side*/, PHG4HitContainer::ConstIterator /*hiter*/, TNtuple * /*ntpad*/, TNtuple * /*nthit*/)=0;// { return {}; }
  void Detector(const std::string &name) { detector = name; }

 protected:
  std::string detector;
};

#endif
