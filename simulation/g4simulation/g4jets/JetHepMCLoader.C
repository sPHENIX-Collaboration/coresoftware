// $Id: $

/*!
 * \file JetHepMCLoader.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "JetHepMCLoader.h"

#include "Jet.h"
#include "JetMapV1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/getClass.h>

#include <TH2F.h>

using namespace std;

JetHepMCLoader::JetHepMCLoader(const std::string &jetInputCategory)
  : SubsysReco("JetHepMCLoader_" + jetInputCategory)
  , m_jetInputCategory(jetInputCategory)
  , m_saveQAPlots(false)
{
}

JetHepMCLoader::~JetHepMCLoader()
{
}

int JetHepMCLoader::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetHepMCLoader::InitRun(PHCompositeNode *topNode)
{
  return CreateNodes(topNode);
}

int JetHepMCLoader::process_event(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetHepMCLoader::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void JetHepMCLoader::addJet(
    const std::string &name,
    int embeddingID,
    const std::string &algorithmName,
    double parameter,
    int tagPID,
    int tagStatus)
{
  hepmc_jet_src src{name, embeddingID, algorithmName, parameter, tagPID, tagStatus};

  m_jetSrc.insert(src);
}

int JetHepMCLoader::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  for (const hepmc_jet_src &src : m_jetSrc)
  {
    // Create the AntiKt node if required
    PHCompositeNode *AlgoNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", src.m_algorithmName.c_str()));
    if (!AlgoNode)
    {
      AlgoNode = new PHCompositeNode(src.m_algorithmName.c_str());
      dstNode->addNode(AlgoNode);
    }

    // Create the Input node if required
    PHCompositeNode *InputNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", m_jetInputCategory.c_str()));
    if (!InputNode)
    {
      InputNode = new PHCompositeNode(m_jetInputCategory.c_str());
      AlgoNode->addNode(InputNode);
    }

    JetMap *jets = findNode::getClass<JetMap>(topNode, src.m_name);
    if (!jets)
    {
      jets = new JetMapV1();
      PHIODataNode<PHObject> *JetMapNode = new PHIODataNode<PHObject>(jets, src.m_name.c_str(), "PHObject");
      InputNode->addNode(JetMapNode);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
