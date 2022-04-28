//CosmicSpray class
// Author: Daniel Lis
// Brief: Particel generator Class that sources a muon with a vertex and momentum that should mimic real life
// modified by Shuhang on 03/2022: now this class serves as a wrapper class that drives "EcoMug" 
#include "CosmicSpray.h"
#include "EcoMug.h"
#include <g4main/PHG4InEvent.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4Utils.h>
#include <g4main/PHG4ParticleGeneratorBase.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>      // for PHDataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <iostream>
#include <TSystem.h>
#include "TROOT.h"
#include "TF3.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>  // for operator<<, endl, basic_ostream
#include <memory>    // for allocator_traits<>::value_type

#include <vector>  // for vector, vector<>::const_iterator

// Declarations
// fixes y plane of cosmic spray here
double CosmicSpray::_y_fix;
// max and min for x spray plane geoemtry
double CosmicSpray::_x_max;
double CosmicSpray::_x_min;
// max and min for z spray plane geometry
double CosmicSpray::_z_max;
double CosmicSpray::_z_min;
// gun energy
double CosmicSpray::_gun_e;





class PHCompositeNode;
class PHG4Particle;
class PHG4ParticleGeneratorBase;


bool CosmicSpray::InDetector(double x, double y, double z){
  double gap = 5;
  if(x > _x_max) return false;
  if(x < -_x_max) return false;
  if(z > _z_max+ gap) return false;
  if(z < -_z_max - gap) return false;
  if(y > _y_fix + gap)return false;
  if(y < -_y_fix - gap) return false;
  return true;
}


CosmicSpray::CosmicSpray(const std::string &name = "COSMICS", const double R = 650, const int &debug = 0)
  : PHG4ParticleGeneratorBase(name)
{
  _x_max = 264.71;
  _x_min = 183.3;
  _z_max = 304.91;
  _z_min = -304.91;
  _y_fix = _x_max;

 
 
  gen.SetUseHSphere(); // half-spherical surface generation
  gen.SetHSphereRadius(R/100); // half-sphere radius
  gen.SetHSphereCenterPosition({{0., 0., -_y_fix / 100}});
  gen.SetMinimumMomentum(0.5);
  _debug = debug;
  _R = R;
  return;
}



int CosmicSpray::process_event(PHCompositeNode *topNode)
{
  // set_vertex
  if(_debug) std::cout<<"Processing Event"<<std::endl;
  std::string pdgname = "mu-";
  int pdgcode  = 13;
  int trackid = 0;
  double gun_t = 0.0;
  double gun_x =0, gun_y =0, gun_z = 0;
  double gun_px = 0, gun_py = 0, gun_pz = 0;
  bool GoodEvent = false;
  while(!GoodEvent){
    gen.Generate();
    std::array<double, 3> muon_position = gen.GetGenerationPosition();
    double muon_p = gen.GetGenerationMomentum();
    _gun_e = sqrt(0.105658* 0.105658 + muon_p * muon_p);
    double tr = gen.GetGenerationTheta();
    double pr = gen.GetGenerationPhi();  
    double muon_charge = gen.GetCharge();
    
    if(muon_charge >0){
      pdgcode = -13;
      pdgname = "mu+";
  }
    
    gun_px = muon_p*sin(tr)*sin(pr);
    gun_py = -1*fabs(muon_p*cos(tr));
    gun_pz = muon_p*sin(tr)*cos(pr);
    
    
    
    gun_x = muon_position[1] * 100;
    gun_y = muon_position[2] * 100;
    gun_z = muon_position[0] * 100;
    
    //unit vectors of muon momentum
    double upx = gun_px/muon_p;
    double upy = gun_py/muon_p;
    double upz = gun_pz/muon_p;
    
    // check if the muon goes in the detector
    
    double x1 = gun_x;
    double y1 = gun_y;
    double z1 = gun_z;
  
    double L = 0;
    
    while(y1 > -_y_fix && L < _R){
      L += 10;
      x1 += upx * 10;
      y1 += upy * 10;
      z1 += upz * 10;
    
      if(InDetector(x1,y1,z1)){
	GoodEvent = true;
	break;
      }
      
    }
    
  }
  
  if(_debug) std::cout<<"Momentum: "<<gun_px<<" / "<<gun_py<<" / "<<gun_pz<<std::endl;
  if(_debug)std::cout<<"total mom: "<<_gun_e<<std::endl;
  if(_debug)std::cout<<"Before adding vertex"<<std::endl;
  _InEvent = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");

  int vtxindex = _InEvent->AddVtx(gun_x, gun_y ,gun_z,gun_t);
  if(_debug)std::cout<<"After adding vertex"<<std::endl;

  PHG4Particle *particle = new PHG4Particlev2();
  particle->set_track_id(trackid);
  if(_debug)std::cout<<"track_id: "<<trackid<<std::endl;
  particle->set_vtx_id(vtxindex);
  if(_debug)std::cout<<"vtxindex: "<<vtxindex<<std::endl;
  particle->set_parent_id(0);
  if(_debug)std::cout<<"parent_id "<<std::endl;
  particle->set_name(pdgname);
  if(_debug)std::cout<<"pdgname: "<<pdgname<<std::endl;
  particle->set_pid(pdgcode);
  if(_debug)std::cout<<"pdgcode: "<<pdgcode<<std::endl;
  particle->set_px(gun_px);
  if(_debug)std::cout<<"px"<<std::endl;
  particle->set_py(gun_py);
  if(_debug)std::cout<<"py"<<std::endl;
  particle->set_pz(gun_pz);
  if(_debug)std::cout<<"pz"<<std::endl;
  particle->set_e(_gun_e);
  if(_debug)std::cout<<"ene"<<std::endl;

  _InEvent->AddParticle(vtxindex, particle);
  if(_debug)std::cout<<"left COSMICS"<<std::endl;

  return 0;
}
