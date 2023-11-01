#include "MbdReco.h"
#include "MbdEvent.h"
#include "MbdPmtContainerV1.h"
#include "MbdVertexMapv1.h"
#include "MbdVertexv1.h"
//#include "MbdOutV2.h"
#include "MbdOutV1.h"
#include "MbdGeomV1.h"

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
MbdReco::MbdReco(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
MbdReco::~MbdReco()
{
}

//____________________________________________________________________________..
int MbdReco::Init(PHCompositeNode *)
{
  m_gaussian = std::make_unique<TF1>("gaussian", "gaus", 0, 20);
  m_gaussian->FixParameter(2, m_tres);

  m_mbdevent = new MbdEvent();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MbdReco::InitRun(PHCompositeNode *topNode)
{
  if (createNodes(topNode) == Fun4AllReturnCodes::ABORTEVENT)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_mbdevent->InitRun();

  int ret = getNodes(topNode);
  return ret;
}

//____________________________________________________________________________..
int MbdReco::process_event(PHCompositeNode *topNode)
{

  getNodes(topNode);

  if ( m_event!=nullptr && m_mbdpmts!=nullptr )
  {
    int status = m_mbdevent->SetRawData( m_event, m_mbdpmts );
    if ( status<0 ) return EVENT_OK; // there wasn't good data in BBC/MBD
  }
  m_mbdevent->Calculate( m_mbdpmts, m_mbdout );

  // For multiple global vertex
  if (m_mbdevent->get_bbcn(0) > 0 && m_mbdevent->get_bbcn(1) > 0)
  {
    auto vertex = std::make_unique<MbdVertexv1>();
    vertex->set_t( m_mbdevent->get_bbct0() );
    vertex->set_z( m_mbdevent->get_bbcz() );
    vertex->set_z_err( 0.6 );
    vertex->set_t_err( m_tres );

    /*
    for (int iarm = 0; iarm < 2; iarm++)
    {
      vertex->set_bbc_ns( iarm, m_mbdevent->get_bbcn(iarm), m_mbdevent->get_bbcq(iarm), m_mbdevent->get_bbct(iarm) );
    }
    */

    m_mbdvtxmap->insert(vertex.release());
  }

  if (Verbosity() > 0)
  {
    std::cout << "mbd vertex z and t0 " << m_mbdevent->get_bbcz() << ", " << m_mbdevent->get_bbct0() << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MbdReco::End(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdReco::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  PHCompositeNode *bbcNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "MBD"));
  if (!bbcNode)
  {
    bbcNode = new PHCompositeNode("MBD");
    dstNode->addNode(bbcNode);
  }

  m_mbdout = findNode::getClass<MbdOut>(bbcNode, "MbdOut");
  if (!m_mbdout)
  {
    m_mbdout = new MbdOutV1();
    PHIODataNode<PHObject> *MbdOutNode = new PHIODataNode<PHObject>(m_mbdout, "MbdOut", "PHObject");
    bbcNode->addNode(MbdOutNode);
  }

  m_mbdpmts = findNode::getClass<MbdPmtContainerV1>(bbcNode, "MbdPmtContainer");
  if (!m_mbdpmts)
  {
    m_mbdpmts = new MbdPmtContainerV1();
    PHIODataNode<PHObject> *MbdPmtContainerNode = new PHIODataNode<PHObject>(m_mbdpmts, "MbdPmtContainer", "PHObject");
    bbcNode->addNode(MbdPmtContainerNode);
  }

  m_mbdvtxmap = findNode::getClass<MbdVertexMap>(bbcNode, "MbdVertexMap");
  if (!m_mbdvtxmap)
  {
    m_mbdvtxmap = new MbdVertexMapv1();
    PHIODataNode<PHObject> *VertexMapNode = new PHIODataNode<PHObject>(m_mbdvtxmap, "MbdVertexMap", "PHObject");
    bbcNode->addNode(VertexMapNode);
  }

  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << PHWHERE << "RUN Node missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  PHNodeIterator runiter(runNode);

  PHCompositeNode *bbcRunNode = dynamic_cast<PHCompositeNode *>(runiter.findFirst("PHCompositeNode", "MBD"));
  if (!bbcRunNode)
  {
    bbcRunNode = new PHCompositeNode("MBD");
    runNode->addNode(bbcRunNode);
  }

  m_mbdgeom = findNode::getClass<MbdGeom>(bbcRunNode, "MbdGeom");
  if (!m_mbdgeom)
  {
    m_mbdgeom = new MbdGeomV1();
    PHIODataNode<PHObject> *MbdGeomNode = new PHIODataNode<PHObject>(m_mbdgeom, "MbdGeom", "PHObject");
    bbcRunNode->addNode(MbdGeomNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdReco::getNodes(PHCompositeNode *topNode)
{
  // Get the bbc prdf data to mpcRawContent
  m_event = findNode::getClass<Event>(topNode,"PRDF");
  //cout << "event addr " << (unsigned int)m_event << endl;

  if ( m_event==nullptr )
  {
    static int counter = 0;
    if ( counter < 1 )
    {
      cout << PHWHERE << "Unable to get PRDF, assuming this is simulation" << endl;
      counter++;
    }
  }

  // MbdPmtContainer
  m_mbdpmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  if (!m_mbdpmts)
  {
    std::cout << PHWHERE << " MbdPmtContainer node not found on node tree" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_mbdvtxmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
  if (!m_mbdvtxmap)
  {
    std::cout << PHWHERE << "MbdVertexMap node not found on node tree" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
