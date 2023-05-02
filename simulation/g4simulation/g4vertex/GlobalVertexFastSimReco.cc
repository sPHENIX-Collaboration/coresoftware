#include "GlobalVertexFastSimReco.h"

#include <globalvertex/GlobalVertex.h>     // for GlobalVertex
#include <globalvertex/GlobalVertexMap.h>  // for GlobalVertexMap
#include <globalvertex/GlobalVertexMapv1.h>
#include <globalvertex/GlobalVertexv1.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_alloc, gsl_rng_free

#include <cmath>
#include <cstdlib>  // for exit
#include <iostream>

GlobalVertexFastSimReco::GlobalVertexFastSimReco(const std::string &name)
  : SubsysReco(name)
{
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
}

GlobalVertexFastSimReco::~GlobalVertexFastSimReco()
{
  gsl_rng_free(RandomGenerator);
}

int GlobalVertexFastSimReco::InitRun(PHCompositeNode *topNode)
{
  if (isnan(_x_smear) ||
      isnan(_y_smear) ||
      isnan(_z_smear) ||
      isnan(_t_smear))
  {
    std::cout << PHWHERE << "::ERROR - smearing must be defined for (x,y,z,t) via set_?_smearing(float)" << std::endl;
    exit(-1);
  }

  unsigned int seed = PHRandomSeed();  // fixed seed handled in PHRandomSeed()
  gsl_rng_set(RandomGenerator, seed);

  if (Verbosity() > 0)
  {
    std::cout << "=================== GlobalVertexFastSimReco::InitRun() ====================" << std::endl;
    std::cout << " x smearing: " << _x_smear << " cm " << std::endl;
    std::cout << " y smearing: " << _y_smear << " cm " << std::endl;
    std::cout << " z smearing: " << _z_smear << " cm " << std::endl;
    std::cout << " t smearing: " << _t_smear << " cm " << std::endl;
    std::cout << " random seed: " << seed << std::endl;
    std::cout << "===========================================================================" << std::endl;
  }

  return CreateNodes(topNode);
}

int GlobalVertexFastSimReco::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "GlobalVertexFastSimReco::process_event -- entered" << std::endl;
  }

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthinfo)
  {
    std::cout << PHWHERE << "::ERROR - cannot find G4TruthInfo" << std::endl;
    exit(-1);
  }

  GlobalVertexMap *vertexes = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexes)
  {
    std::cout << PHWHERE << "::ERROR - cannot find GlobalVertexMap" << std::endl;
    exit(-1);
  }

  //---------------------
  // Smear Truth Vertexes (only one per crossing right now)
  //---------------------

  PHG4VtxPoint *point = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());

  GlobalVertex *vertex = new GlobalVertexv1(GlobalVertex::TRUTH);
  vertex->set_x(point->get_x());
  vertex->set_y(point->get_y());
  vertex->set_z(point->get_z());
  vertex->set_t(point->get_t());
  vertex->set_t_err(0.);
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      vertex->set_error(i, j, 0.0);
    }
  }
  vertexes->insert(vertex);

  vertex = new GlobalVertexv1(GlobalVertex::SMEARED);

  vertex->set_x(point->get_x() + gsl_ran_gaussian(RandomGenerator, _x_smear));
  vertex->set_y(point->get_y() + gsl_ran_gaussian(RandomGenerator, _y_smear));
  vertex->set_z(point->get_z() + gsl_ran_gaussian(RandomGenerator, _z_smear));

  vertex->set_error(0, 0, _x_smear * _x_smear);
  vertex->set_error(0, 1, 0.0);
  vertex->set_error(0, 2, 0.0);

  vertex->set_error(1, 0, 0.0);
  vertex->set_error(1, 1, _y_smear * _y_smear);
  vertex->set_error(1, 2, 0.0);

  vertex->set_error(2, 0, 0.0);
  vertex->set_error(2, 1, 0.0);
  vertex->set_error(2, 2, _z_smear * _z_smear);

  vertex->set_t(point->get_t() + gsl_ran_gaussian(RandomGenerator, _t_smear));
  vertex->set_t_err(_t_smear);

  vertexes->insert(vertex);

  return Fun4AllReturnCodes::EVENT_OK;
}

int GlobalVertexFastSimReco::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // store the GLOBAL stuff under a sub-node directory
  PHCompositeNode *globalNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!globalNode)
  {
    globalNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(globalNode);
  }

  // create the GlobalVertexMap
  GlobalVertexMap *vertexes = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexes)
  {
    vertexes = new GlobalVertexMapv1();
    PHIODataNode<PHObject> *VertexMapNode = new PHIODataNode<PHObject>(vertexes, "GlobalVertexMap", "PHObject");
    globalNode->addNode(VertexMapNode);
  }
  else
  {
    std::cout << PHWHERE << "::ERROR - GlobalVertexMap pre-exists, but should not if running FastSim" << std::endl;
    exit(-1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
