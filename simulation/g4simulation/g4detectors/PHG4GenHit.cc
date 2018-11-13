#include "PHG4GenHit.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <cmath>

using namespace std;

PHG4GenHit::PHG4GenHit(const string &name):
  SubsysReco(name),
  phi(NAN),
  theta(NAN),
  eloss(NAN),
  layer(-9999)
{}

int
PHG4GenHit::process_event(PHCompositeNode *topNode)
{
  string hitnodename = "G4HIT_" + detector;
  string geonodename = "CYLINDERGEOM_" + detector;
  PHG4CylinderGeomContainer *geo =  findNode::getClass<PHG4CylinderGeomContainer>(topNode , geonodename.c_str());
  if (! geo)
    {
      cout << "cannot find geo node " << geonodename << endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  PHG4HitContainer *hits_ = findNode::getClass<PHG4HitContainer>( topNode , hitnodename.c_str());
  if (! hits_)
    {
      cout << "cannot find hit node " << hitnodename << endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  PHG4CylinderGeom *mygeom = geo->GetLayerGeom(layer);
  double inner_radius = mygeom->get_radius();
  double outer_radius = inner_radius + mygeom->get_thickness();
  PHG4Hit *hit = new PHG4Hitv1();
  hit->set_layer((unsigned int)layer);
  double x0 = inner_radius * cos(phi * M_PI / 180.);
  double y0 = inner_radius * sin(phi * M_PI / 180.);
  double z0 = inner_radius * cos(theta * M_PI / 180.);
  double x1 = outer_radius * cos(phi * M_PI / 180.);
  double y1 = outer_radius * sin(phi * M_PI / 180.);
  double z1 = outer_radius * cos(theta * M_PI / 180.);
  hit->set_x(0, x0);
  hit->set_y(0, y0);
  hit->set_z(0, z0);
  hit->set_x(1, x1);
  hit->set_y(1, y1);
  hit->set_z(1, z1);
  hit->set_edep(eloss);
  hit->set_trkid(-1);
  hits_->AddHit(layer, hit);
  if (Verbosity() > 0)
    {
  cout << "phi " << phi << " inner rad: " << inner_radius
       << ", outer rad: " << outer_radius
       << " x0/y0/z0: " << x0 << "/" << y0 << "/" << z0
       << " x1/y1/z1: " << x1 << "/" << y1 << "/" << z1
         << " edep: " << eloss
         << endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
}
