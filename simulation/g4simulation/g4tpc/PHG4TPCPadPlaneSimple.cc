#include "PHG4TPCPadPlaneSimple.h"
#include "PHG4CellTPCv1.h"

#include <g4detectors/PHG4CellContainer.h>

#include <cmath>
#include <iostream>

using namespace std;

PHG4TPCPadPlaneSimple::PHG4TPCPadPlaneSimple(const string &name):
PHG4TPCPadPlane(name),
max_active_radius(NAN),
min_active_radius(NAN),
rbinwidth(NAN),
phibinwidth(NAN),
tbinwidth(NAN)
{
  InitializeParameters();
  return;
}


void PHG4TPCPadPlaneSimple::MapToPadPlane(PHG4CellContainer *g4cells, const double x_gem, const double y_gem, const double t_gem,  PHG4HitContainer::ConstIterator hiter)
{
  double phi = atan2(y_gem,x_gem);
  double rad_gem = sqrt(x_gem*x_gem + y_gem*y_gem);
  int phibin = (phi+M_PI)/phibinwidth;
  int radbin = (rad_gem-min_active_radius)/rbinwidth;
  int tbin = t_gem/tbinwidth;
//  ntpad->Fill(t_gem,phi,rad_gem,phibin,radbin);
  // cout << ", phi: " << phi
  //      << ", phibin: " << phibin
  //      << ", rad: " << rad_gem
  //      << ", radbin: " << radbin
  //      << ", t_gem: " << t_gem 
  //      << endl;
  PHG4CellDefs::keytype key = PHG4CellDefs::TPCBinning::genkey(0,radbin,phibin);
  PHG4Cell *cell = g4cells->findCell(key);
  if (! cell)
  {
    cell = new PHG4CellTPCv1(key);
    g4cells->AddCell(cell);
  }
  cell->add_edep(key,tbin,1.);
   return;
}

void PHG4TPCPadPlaneSimple::SetDefaultParameters()
{
  set_default_int_param("nbins_phi",360);
  set_default_int_param("nbins_r",40);

  set_default_double_param("binwidth_t",53.); //  2*Rhic clock = 106/2.
  set_default_double_param("max_active_radius",75.);
  set_default_double_param("min_active_radius",30.);
  return;
}

void PHG4TPCPadPlaneSimple::UpdateInternalParameters()
{
  max_active_radius = get_double_param("max_active_radius");
  min_active_radius = get_double_param("min_active_radius");
  rbinwidth = (max_active_radius - min_active_radius)/get_int_param("nbins_r");
  phibinwidth = 2*M_PI/get_int_param("nbins_phi");
  tbinwidth = get_double_param("binwidth_t");
}
