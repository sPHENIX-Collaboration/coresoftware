#ifndef PHG4TPCPadPlaneSimple_h
#define PHG4TPCPadPlaneSimple_h

#include "PHG4TPCPadPlane.h"

class PHG4CellContainer;

class PHG4TPCPadPlaneSimple: public PHG4TPCPadPlane
{
public:
  PHG4TPCPadPlaneSimple(const std::string& name = "SimplePadPlane");
  virtual ~PHG4TPCPadPlaneSimple(){}

  void MapToPadPlane(PHG4CellContainer *g4cells, const double x_gem, const double y_gem, const double t_gem, PHG4HitContainer::ConstIterator hiter);

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
