#include "PHG4InEventReadBack.h"
#include "PHG4InEvent.h"
#include "PHG4Particle.h"                 // for PHG4Particle
#include "PHG4Particlev1.h"
#include "PHG4VtxPointv1.h"

#include <vararray/VariableArray.h>
#include <vararray/VariableArrayUtils.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                 // for PHNode
#include <phool/PHNodeIterator.h>         // for PHNodeIterator
#include <phool/PHObject.h>               // for PHObject

#include <iostream>

using namespace std;

PHG4InEventReadBack::PHG4InEventReadBack(const std::string &name): SubsysReco(name)
{}

int
PHG4InEventReadBack::InitRun(PHCompositeNode *topNode)
{
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode,"PHG4INEVENT");
  if (!ineve)
    {
      PHNodeIterator iter( topNode );
      PHCompositeNode *dstNode;
      dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

      ineve = new PHG4InEvent();
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(ineve, "PHG4INEVENT", "PHObject");
      dstNode->addNode(newNode);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int
PHG4InEventReadBack::process_event(PHCompositeNode *topNode)
{
  PHG4InEvent *inEvent = findNode::getClass<PHG4InEvent>(topNode,"PHG4INEVENT");
  if (!inEvent)
    {
      cout << "no PHG4INEVENT node found" << endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  VariableArray *vtxarray = findNode::getClass<VariableArray>(topNode,"PHG4Vtx_VarArray");
  if (!vtxarray)
    {
      cout << "no PHG4Vtx_VarArray node found" << endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }

  VariableArray *particlearray = findNode::getClass<VariableArray>(topNode,"PHG4Particle_VarArray");
  if (!particlearray)
    {
      cout << "no PHG4Particle_VarArray node found" << endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  inEvent->Reset();
  unsigned int size = vtxarray->get_array_size();
  const short int *sval = vtxarray->get_array();
  PHG4VtxPointv1 vtx;
  while(size > 0)
    {
      int vtxid = *sval++;
      size --;
      vtx.set_x(VariableArrayUtils::ShortBitsToFloat(*sval++));
      size --;
      vtx.set_y(VariableArrayUtils::ShortBitsToFloat(*sval++));
      size --;
      vtx.set_z(VariableArrayUtils::ShortBitsToFloat(*sval++));
      size--;
      vtx.set_t(VariableArrayUtils::ShortBitsToFloat(*sval++));
      size--;
      inEvent->AddVtx(vtxid,vtx);
     }

  size = particlearray->get_array_size();
  sval = particlearray->get_array();
  while(size > 0)
    {
      PHG4Particle *particle = new PHG4Particlev1();
      int vtxid = *sval++;
      size --;
      particle->set_pid(*sval++);
      size --;
      particle->set_px(VariableArrayUtils::ShortBitsToFloat(*sval++));
      size --;
      particle->set_py(VariableArrayUtils::ShortBitsToFloat(*sval++));
      size --;
      particle->set_pz(VariableArrayUtils::ShortBitsToFloat(*sval++));
      size --;
      inEvent->AddParticle(vtxid, particle);
    }
  //  inEvent->identify();
  return Fun4AllReturnCodes::EVENT_OK;
}

int
PHG4InEventReadBack::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
