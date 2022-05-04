//CosmicSpray class
// Author: Daniel Lis
// Brief: Particel generator Class that sources a muon with a vertex and momentum that should mimic real life
// modified by Shuhang on 03/2022: now this class serves as a wrapper class that drives "EcoMug"

//For the muon rate calculation, if using the default setting: time(second) = nevents * 1.434e-4, then scale the histogram by 1/time to get the rate.

#include "CosmicSpray.h"
#include "EcoMug.h"

#include "PHG4InEvent.h"
#include "PHG4Particle.h"
#include "PHG4Particlev2.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>      // for PHDataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/PHRandomSeed.h>


#include <array>  // for array
#include <cmath>
#include <iostream>  // for operator<<, endl, basic_ostream

bool CosmicSpray::InDetector(double x, double y, double z)
{
  double gap = 5;
  if (x > _x_max) return false;
  if (x < -_x_max) return false;
  if (z > _z_max + gap) return false;
  if (z < -_z_max - gap) return false;
  if (y > _y_fix + gap) return false;
  if (y < -_y_fix - gap) return false;
  return true;
}

CosmicSpray::CosmicSpray(const std::string &name, const double R)
  : SubsysReco(name)
{
  _x_max = 264.71;
  _x_min = 183.3;
  _z_max = 304.91;
  _z_min = -304.91;
  _y_fix = _x_max;

  unsigned int seed = PHRandomSeed();
  gen.SetSeed(seed);
  gen.SetUseHSphere();            // half-spherical surface generation
  gen.SetHSphereRadius(R / 100);  // half-sphere radius
  gen.SetHSphereCenterPosition({{0., 0., -_y_fix / 100}});
  gen.SetMinimumMomentum(0.5);
  _R = R;
  return;
}

int CosmicSpray::InitRun(PHCompositeNode *topNode)
{
  PHG4InEvent *inevent = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!inevent)
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode;
    dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

    inevent = new PHG4InEvent();
    PHDataNode<PHObject> *newNode = new PHDataNode<PHObject>(inevent, "PHG4INEVENT", "PHObject");
    dstNode->addNode(newNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int CosmicSpray::process_event(PHCompositeNode *topNode)
{
  // set_vertex
  if (Verbosity() > 0) std::cout << "Processing Event" << std::endl;
  std::string pdgname = "mu-";
  int pdgcode = 13;
  int trackid = 0;
  double gun_t = 0.0;
  double gun_x = 0, gun_y = 0, gun_z = 0;
  double gun_px = 0, gun_py = 0, gun_pz = 0;
  bool GoodEvent = false;
  while (!GoodEvent)
  {
    gen.Generate();
    std::array<double, 3> muon_position = gen.GetGenerationPosition();
    double muon_p = gen.GetGenerationMomentum();
    _gun_e = sqrt(0.105658 * 0.105658 + muon_p * muon_p);
    double tr = gen.GetGenerationTheta();
    double pr = gen.GetGenerationPhi();
    double muon_charge = gen.GetCharge();

    if (muon_charge > 0)
    {
      pdgcode = -13;
      pdgname = "mu+";
    }

    gun_px = muon_p * sin(tr) * sin(pr);
    gun_py = -1 * fabs(muon_p * cos(tr));
    gun_pz = muon_p * sin(tr) * cos(pr);

    gun_x = muon_position[1] * 100;
    gun_y = muon_position[2] * 100;
    gun_z = muon_position[0] * 100;

    //unit vectors of muon momentum
    double upx = gun_px / muon_p;
    double upy = gun_py / muon_p;
    double upz = gun_pz / muon_p;

    // check if the muon goes in the detector

    double x1 = gun_x;
    double y1 = gun_y;
    double z1 = gun_z;

    double L = 0;
    

    while (y1 > -_y_fix && L < _R)
    {
      L += 10;
      x1 += upx * 10;
      y1 += upy * 10;
      z1 += upz * 10;

      if (InDetector(x1, y1, z1))
      {
        GoodEvent = true;
        break;
      }
    }
  }

  PHG4InEvent *inevent = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  int vtxindex = inevent->AddVtx(gun_x, gun_y, gun_z, gun_t);
  
  PHG4Particle *particle = new PHG4Particlev2();
  particle->set_track_id(trackid);
  particle->set_vtx_id(vtxindex);
  particle->set_parent_id(0);
  particle->set_name(pdgname);
  particle->set_pid(pdgcode);
  particle->set_px(gun_px);
  particle->set_py(gun_py);
  particle->set_pz(gun_pz);
  particle->set_e(_gun_e);
  inevent->AddParticle(vtxindex, particle);
  if (Verbosity() > 0)
  {
    std::cout << "Momentum: " << gun_px << " / " << gun_py << " / " << gun_pz << std::endl;
    std::cout << "total mom: " << _gun_e << std::endl;
    std::cout << "track_id: " << trackid << std::endl;
    std::cout << "vtxindex: " << vtxindex << std::endl;
    std::cout << "pdgname: " << pdgname << std::endl;
    std::cout << "pdgcode: " << pdgcode << std::endl;
  }
  return 0;
}
