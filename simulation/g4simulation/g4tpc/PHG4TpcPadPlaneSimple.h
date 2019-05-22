#ifndef G4TPC_PHG4TPCPADPLANESIMPLE_H
#define G4TPC_PHG4TPCPADPLANESIMPLE_H

#include "PHG4TpcPadPlane.h"

class PHG4CellContainer;

class PHG4TpcPadPlaneSimple : public PHG4TpcPadPlane
{
 public:
  PHG4TpcPadPlaneSimple(const std::string& name = "PHG4TpcPadPlaneSimple");
  virtual ~PHG4TpcPadPlaneSimple() {}

  void MapToPadPlane(PHG4CellContainer* g4cells, const double x_gem, const double y_gem, const double t_gem, PHG4HitContainer::ConstIterator hiter);

  void SetDefaultParameters();
  void UpdateInternalParameters();

 protected:
  double max_active_radius;
  double min_active_radius;
  double rbinwidth;
  double phibinwidth;
  double tbinwidth;
};

#endif
