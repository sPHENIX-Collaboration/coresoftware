#include "PHG4InEventCompress.h"

#include "PHG4InEvent.h"
#include "PHG4VtxPoint.h"
#include "PHG4Particle.h"

#include <vararray/VariableArray.h>
#include <vararray/VariableArrayIds.h>
#include <vararray/VariableArrayUtils.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                 // for PHNode
#include <phool/PHNodeIterator.h>         // for PHNodeIterator
#include <phool/PHObject.h>               // for PHObject

#include <cstdlib>
#include <iostream>
#include <map>                            // for _Rb_tree_const_iterator
#include <utility>                        // for pair
#include <vector>

using namespace std;

PHG4InEventCompress::PHG4InEventCompress(const std::string &name): 
  SubsysReco(name),
  vtxarray(nullptr),
  particlearray(nullptr)
{}

int
PHG4InEventCompress::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));


  PHIODataNode<PHObject> *PHObjectIONode;

  vtxarray = new VariableArray(varids::G4VTXV1);
  PHObjectIONode = new PHIODataNode<PHObject>(vtxarray, "PHG4Vtx_VarArray", "PHObject");
  dstNode->addNode(PHObjectIONode);

  particlearray = new VariableArray(varids::G4PARTICLEV1);
  PHObjectIONode = new PHIODataNode<PHObject>(particlearray, "PHG4Particle_VarArray", "PHObject");
  dstNode->addNode(PHObjectIONode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int
PHG4InEventCompress::process_event(PHCompositeNode *topNode)
{
  PHG4InEvent *inEvent = findNode::getClass<PHG4InEvent>(topNode,"PHG4INEVENT");
  if (!inEvent)
    {
      cout << "no PHG4INEVENT node found" << endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }

  map<int, PHG4VtxPoint *>::const_iterator vtxiter;
  std::pair< std::map<int, PHG4VtxPoint *>::const_iterator, std::map<int, PHG4VtxPoint *>::const_iterator > vtxbegin_end = inEvent->GetVertices();
  vector<short> svtxvec;
  for (vtxiter = vtxbegin_end.first; vtxiter != vtxbegin_end.second; ++vtxiter)
    {
      if ((*vtxiter).first > 0xFFFF)
	{
	  cout << "id of vertex " << (*vtxiter).first << " exceeds max val of " << 0xFFFF << endl;
	  exit(1);
	}
      svtxvec.push_back((*vtxiter).first);
      svtxvec.push_back(VariableArrayUtils::FloatToShortBits((*vtxiter->second).get_x()));
      svtxvec.push_back(VariableArrayUtils::FloatToShortBits((*vtxiter->second).get_y()));
      svtxvec.push_back(VariableArrayUtils::FloatToShortBits((*vtxiter->second).get_z()));
      svtxvec.push_back(VariableArrayUtils::FloatToShortBits((*vtxiter->second).get_t()));
    }
  vtxarray->set_val(svtxvec);

  pair<multimap<int, PHG4Particle *>::const_iterator, multimap<int, PHG4Particle *>::const_iterator > particlebegin_end = inEvent->GetParticles();
  multimap<int,PHG4Particle *>::const_iterator particle_iter;
  vector<short> spartvec;
  for (particle_iter = particlebegin_end.first; particle_iter != particlebegin_end.second; ++particle_iter)
    {
      if ((*particle_iter).first > 0xFFFF)
	{
	  cout << "id of vertex " << (*particle_iter).first << " exceeds max val of " << 0xFFFF << endl;
	  exit(1);
	}
      spartvec.push_back((*particle_iter).first);
      if (abs((*particle_iter->second).get_pid()) > 0xFFFF)
	{
	  cout << "pdg code of particle " << (*particle_iter->second).get_pid() << " exceeds max val of " << 0xFFFF << endl;
	  exit(1);
	}
      spartvec.push_back((*particle_iter->second).get_pid());

      spartvec.push_back(VariableArrayUtils::FloatToShortBits((*particle_iter->second).get_px()));
      spartvec.push_back(VariableArrayUtils::FloatToShortBits((*particle_iter->second).get_py()));
      spartvec.push_back(VariableArrayUtils::FloatToShortBits((*particle_iter->second).get_pz()));
    }
  particlearray->set_val(spartvec);
  //  inEvent->identify();
  return Fun4AllReturnCodes::EVENT_OK;
}

int
PHG4InEventCompress::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
