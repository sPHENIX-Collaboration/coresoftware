#include "CaloVtxReco.h"

#include <globalvertex/CaloVertexMapv1.h>
#include <globalvertex/CaloVertexv1.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include "TMath.h"
#include <cmath>
/*
float radius_EM = 93.5;
float radius_OH = 225.87;
*/
//____________________________________________________________________________..
CaloVtxReco::CaloVtxReco(const std::string &name)
  : SubsysReco(name)
{
}

int CaloVtxReco::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << PHWHERE << "RUN Node missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHNodeIterator dstiter(dstNode);

  PHCompositeNode *globalNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!globalNode)
  {
    globalNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(globalNode);
  }

  m_calovtxmap = findNode::getClass<CaloVertexMap>(globalNode, "CaloVertexMap");
  if (!m_calovtxmap)
  {
    m_calovtxmap = new CaloVertexMapv1();
    PHIODataNode<PHObject> *VertexMapNode = new PHIODataNode<PHObject>(m_calovtxmap, "CaloVertexMap", "PHObject");
    globalNode->addNode(VertexMapNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloVtxReco::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "Initializing!" << std::endl;
  }
  int status = createNodes(topNode);
  if (status)
  {
    return status;
  }

  for (auto &algo : m_algos)
    {
      status = algo->Init(topNode);
      if (status)
	{
	  return status;
	}
    }
   
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloVtxReco::process_event(PHCompositeNode *topNode)
{
  

  for (auto &algo : m_algos)
    {
      float tempz = std::numeric_limits<float>::quiet_NaN();
      algo->CalculateVertex(topNode, tempz);

      CaloVertex *vertex = new CaloVertexv1();  

      vertex->set_z(tempz);
      vertex->set_z_err(0);
      vertex->set_t(0);
      vertex->set_t_err(0);
      vertex->set_calo_algo(algo->Algo());
      m_calovtxmap->insert(vertex);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}
