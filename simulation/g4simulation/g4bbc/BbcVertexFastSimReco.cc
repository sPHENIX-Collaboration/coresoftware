
#include "BbcVertexFastSimReco.h"

#include "BbcVertexMap.h"                   // for BbcVertexMap
#include "BbcVertexMapv1.h"
#include "BbcVertex.h"                      // for BbcVertex
#include "BbcVertexv1.h"

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>             // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                   // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                 // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>                    // for PHWHERE

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>                    // for gsl_rng_uniform_pos, gsl_...

#include <cstdlib>                         // for exit
#include <cmath>
#include <iostream>

using namespace std;

BbcVertexFastSimReco::BbcVertexFastSimReco(const string &name)
  : SubsysReco(name)
  , m_T_Smear(NAN)
  , m_Z_Smear(NAN)
{
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
}

BbcVertexFastSimReco::~BbcVertexFastSimReco()
{
  gsl_rng_free(RandomGenerator);
}

int BbcVertexFastSimReco::Init(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int BbcVertexFastSimReco::InitRun(PHCompositeNode *topNode)
{
  if (isnan(m_T_Smear) || isnan(m_Z_Smear))
  {
    cout << PHWHERE << "::ERROR - smearing must be defined via setm_T_Smearing(float) and set_z_smeaering(float)" << endl;
    exit(-1);
  }

  unsigned int seed = PHRandomSeed();  // fixed random seed handled in PHRandomSeed()
  gsl_rng_set(RandomGenerator, seed);

  if (Verbosity() > 0)
  {
    cout << "===================== BbcVertexFastSimReco::InitRun() =====================" << endl;
    cout << " t smearing: " << m_T_Smear << " cm " << endl;
    cout << "  z smearing: " << m_Z_Smear << " cm " << endl;
    cout << " random seed: " << seed << endl;
    cout << "===========================================================================" << endl;
  }

  return CreateNodes(topNode);
}

int BbcVertexFastSimReco::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 1) cout << "BbcVertexFastSimReco::process_event -- entered" << endl;

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthinfo)
  {
    cout << PHWHERE << "::ERROR - cannot find G4TruthInfo" << endl;
    exit(-1);
  }

  BbcVertexMap *vertexes = findNode::getClass<BbcVertexMap>(topNode, "BbcVertexMap");
  if (!vertexes)
  {
    cout << PHWHERE << "::ERROR - cannot find BbcVertexMap" << endl;
    exit(-1);
  }

  //---------------------
  // Smear Truth Vertexes (only one per crossing right now)
  //---------------------

  PHG4VtxPoint *point = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
  if (!point) return Fun4AllReturnCodes::EVENT_OK;

  BbcVertex *vertex = new BbcVertexv1();

  if (m_T_Smear >= 0.0)
  {
    vertex->set_t(point->get_t() + gsl_ran_gaussian(RandomGenerator, m_T_Smear));
    vertex->set_t_err(m_T_Smear);
  }
  else
  {
    vertex->set_t(point->get_t() + 2.0 * gsl_rng_uniform_pos(RandomGenerator) * m_T_Smear);
    vertex->set_t_err(fabs(m_T_Smear) / sqrt(12));
  }

  if (m_Z_Smear >= 0.0)
  {
    vertex->set_z(point->get_z() + gsl_ran_gaussian(RandomGenerator, m_Z_Smear));
    vertex->set_z_err(m_Z_Smear);
  }
  else
  {
    vertex->set_z(point->get_z() + 2.0 * gsl_rng_uniform_pos(RandomGenerator) * m_Z_Smear);
    vertex->set_z_err(fabs(m_Z_Smear) / sqrt(12));
  }

  vertexes->insert(vertex);

  return Fun4AllReturnCodes::EVENT_OK;
}

int BbcVertexFastSimReco::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int BbcVertexFastSimReco::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // store the BBC stuff under a sub-node directory
  PHCompositeNode *bbcNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "BBC"));
  if (!bbcNode)
  {
    bbcNode = new PHCompositeNode("BBC");
    dstNode->addNode(bbcNode);
  }

  // create the BbcVertexMap
  BbcVertexMap *vertexes = findNode::getClass<BbcVertexMap>(topNode, "BbcVertexMap");
  if (!vertexes)
  {
    vertexes = new BbcVertexMapv1();
    PHIODataNode<PHObject> *VertexMapNode = new PHIODataNode<PHObject>(vertexes, "BbcVertexMap", "PHObject");
    bbcNode->addNode(VertexMapNode);
  }
  else
  {
    cout << PHWHERE << "::ERROR - BbcVertexMap pre-exists, but should not if running FastSim" << endl;
    exit(-1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
