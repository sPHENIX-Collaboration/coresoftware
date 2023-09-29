#include "BbcReco.h"
#include "BbcDefs.h"
#include "BbcEvent.h"
#include "BbcPmtInfoContainerV1.h"
#include "BbcReturnCodes.h"
#include "BbcVertexMapv1.h"
#include "BbcVertexv2.h"
#include "BbcOutV1.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <Event/Event.h>

#include <TF1.h>
#include <TH1.h>

using namespace std;
using namespace Fun4AllReturnCodes;

//____________________________________________________________________________..
BbcReco::BbcReco(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
BbcReco::~BbcReco()
{
}

//____________________________________________________________________________..
int BbcReco::Init(PHCompositeNode *)
{
  m_gaussian = std::make_unique<TF1>("gaussian", "gaus", 0, 20);
  m_gaussian->FixParameter(2, m_tres);

  m_bbcevent = new BbcEvent();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int BbcReco::InitRun(PHCompositeNode *topNode)
{
  if (createNodes(topNode) == Fun4AllReturnCodes::ABORTEVENT)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  int ret = getNodes(topNode);
  return ret;
}

//____________________________________________________________________________..
int BbcReco::process_event(PHCompositeNode *)
{

  if ( m_event!=0 && m_bbcpmts!=0 ) m_bbcevent->SetRawData( m_event, m_bbcpmts );
  m_bbcevent->Calculate( m_bbcpmts, m_bbcout );

  // This seems to exactly duplicate bbcout. Why?
  if (m_bbcevent->get_bbcn(0) > 0 && m_bbcevent->get_bbcn(1) > 0)
  {
    auto vertex = std::make_unique<BbcVertexv2>();
    vertex->set_t( m_bbcevent->get_bbct0() );
    vertex->set_z( m_bbcevent->get_bbcz() );
    vertex->set_z_err( 0.6 );
    vertex->set_t_err( m_tres );

    for (int iarm = 0; iarm < 2; iarm++)
    {
      vertex->set_bbc_ns( iarm, m_bbcevent->get_bbcn(iarm), m_bbcevent->get_bbcq(iarm), m_bbcevent->get_bbct(iarm) );
    }

    m_bbcvertexmap->insert(vertex.release());
  }

  if (Verbosity() > 0)
  {
    std::cout << "bbc vertex z and t0 " << m_bbcevent->get_bbcz() << ", " << m_bbcevent->get_bbct0() << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int BbcReco::End(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int BbcReco::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  PHCompositeNode *bbcNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "BBC"));
  if (!bbcNode)
  {
    bbcNode = new PHCompositeNode("BBC");
    dstNode->addNode(bbcNode);
  }

  m_bbcout = findNode::getClass<BbcOut>(bbcNode, "BbcOut");
  if (!m_bbcout)
  {
    m_bbcout = new BbcOutV1();
    PHIODataNode<PHObject> *BbcOutNode = new PHIODataNode<PHObject>(m_bbcout, "BbcOut", "PHObject");
    bbcNode->addNode(BbcOutNode);
  }

  m_bbcpmts = findNode::getClass<BbcPmtInfoContainerV1>(topNode, "BbcPmtInfoContainer");
  if (!m_bbcpmts)
  {
    m_bbcpmts = new BbcPmtInfoContainerV1();
    PHIODataNode<PHObject> *BbcPmtInfoContainerNode = new PHIODataNode<PHObject>(m_bbcpmts, "BbcPmtInfoContainer", "PHObject");
    bbcNode->addNode(BbcPmtInfoContainerNode);
  }

  m_bbcvertexmap = findNode::getClass<BbcVertexMap>(bbcNode, "BbcVertexMap");
  if (!m_bbcvertexmap)
  {
    m_bbcvertexmap = new BbcVertexMapv1();
    PHIODataNode<PHObject> *VertexMapNode = new PHIODataNode<PHObject>(m_bbcvertexmap, "BbcVertexMap", "PHObject");
    bbcNode->addNode(VertexMapNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int BbcReco::getNodes(PHCompositeNode *topNode)
{
  // Get the bbc prdf data to mpcRawContent
  m_event = findNode::getClass<Event>(topNode,"PRDF");
  //cout << "event addr " << (unsigned int)m_event << endl;

  if ( m_event==0 )
  {
    static int counter = 0;
    if ( counter < 1 )
    {
      cout << PHWHERE << "Unable to get PRDF, assuming this is simulation" << endl;
      counter++;
    }
  }

  // BbcPmtContainer
  m_bbcpmts = findNode::getClass<BbcPmtInfoContainerV1>(topNode, "BbcPmtInfoContainer");
  if (!m_bbcpmts)
  {
    std::cout << PHWHERE << " BbcPmtInfoContainer node not found on node tree" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_bbcvertexmap = findNode::getClass<BbcVertexMap>(topNode, "BbcVertexMap");
  if (!m_bbcvertexmap)
  {
    std::cout << PHWHERE << "BbcVertexMap node not found on node tree" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
