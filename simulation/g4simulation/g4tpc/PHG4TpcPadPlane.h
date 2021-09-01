#ifndef G4TPC_PHG4TPCPADPLANE_H
#define G4TPC_PHG4TPCPADPLANE_H

#include <fun4all/SubsysReco.h>

#include <g4main/PHG4HitContainer.h>

#include <phparameter/PHParameterInterface.h>

#include <string>                              // for string

class PHG4CellContainer;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class PHCompositeNode;
class PHG4CylinderCellGeomContainer;
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
  virtual int CreateReadoutGeometry(PHCompositeNode *, PHG4CylinderCellGeomContainer *) { return 0; }
  virtual void UpdateInternalParameters() { return; }
  virtual void MapToPadPlane(PHG4CellContainer *, const double, const double, const double, PHG4HitContainer::ConstIterator, TNtuple *, TNtuple *) {}
  virtual void MapToPadPlane(TrkrHitSetContainer *, TrkrHitSetContainer *, TrkrHitTruthAssoc *, const double, const double, const double, PHG4HitContainer::ConstIterator, TNtuple *, TNtuple *) {}
  void Detector(const std::string &name) { detector = name; }

 protected:
  std::string detector;
};

#endif
