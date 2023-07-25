#include "BbcReconstruction.h"
#include "BbcDefs.h"
#include "BbcPmtContainer.h"
#include "BbcReturnCodes.h"
#include "BbcVertexMapv1.h"
#include "BbcVertexv2.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TF1.h>
#include <TH1.h>

//____________________________________________________________________________..
BbcReconstruction::BbcReconstruction(const std::string &name)
  : SubsysReco(name)
{
  h_evt_bbct[0] = nullptr;
  h_evt_bbct[1] = nullptr;
}

//____________________________________________________________________________..
BbcReconstruction::~BbcReconstruction()
{
}

//____________________________________________________________________________..
int BbcReconstruction::Init(PHCompositeNode *)
{
  m_gaussian = std::make_unique<TF1>("gaussian", "gaus", 0, 20);
  m_gaussian->FixParameter(2, m_tres);

  TString name, title;
  for (int iarm = 0; iarm < 2; iarm++)
  {
    //
    name = "hevt_bbct";
    name += iarm;
    title = "bbc times, arm ";
    title += iarm;
    h_evt_bbct[iarm] = new TH1F(name, title, 200, 7.5, 11.5);
    h_evt_bbct[iarm]->SetLineColor(4);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int BbcReconstruction::InitRun(PHCompositeNode *topNode)
{
  if (createNodes(topNode) == Fun4AllReturnCodes::ABORTEVENT)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  int ret = getNodes(topNode);
  return ret;
}

//____________________________________________________________________________..
int BbcReconstruction::process_event(PHCompositeNode *)
{

  std::vector<float> hit_times[2];
  float bbcq[2] = {0};
  float bbcn[2] = {0};
  float bbct[2] = {0};

  for (int ich = 0; ich < 128; ich++)
  {
    const auto pmt = m_bbcpmts->get_pmt(ich);
    const int arm = pmt / 64;
    auto bbchit = m_bbcpmts->get_tdc0(ich);

    if (std::isnan(m_bbcpmts->get_adc(ich)))
    {
      //! PMT with no ADC, skip it
      continue;
    }

    if (bbchit < 9999.)
    {
      //! total charge
      bbcq[arm] += m_bbcpmts->get_adc(ich);

      //! number of hit pmts
      bbcn[arm]++;

      h_evt_bbct[arm]->Fill(bbchit);
      hit_times[arm].push_back(bbchit);
    }
    if (Verbosity() > 0)
    {
      std::cout << "Channel " << ich << " with pmt " << pmt << " in arm " << arm << "  has " << m_bbcpmts->get_adc(ich) << std::endl;
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "nbbc arm 1,2: " << bbcn[0] << ", " << bbcn[1] << std::endl;
  }

  if (bbcn[0] > 0 && bbcn[1] > 0)
  {
    for (int iarm = 0; iarm < 2; iarm++)
    {
      if (hit_times[iarm].empty())
      {
        continue;
      }

      std::sort(hit_times[iarm].begin(), hit_times[iarm].end());
      float earliest = hit_times[iarm][0];
      m_gaussian->SetParameter(0, 5);
      m_gaussian->SetParameter(1, earliest);
      m_gaussian->SetRange(6, earliest + 5 * m_tres);
      h_evt_bbct[iarm]->Fit(m_gaussian.get(), "BLRNQ");
      bbct[iarm] = m_gaussian->GetParameter(1);
    }

    float bbcz = (bbct[0] - bbct[1]) * BbcDefs::C / 2.;
    float bbct0 = (bbct[0] + bbct[1]) / 2.;

    auto vertex = std::make_unique<BbcVertexv2>();
    vertex->set_t(bbct0);
    vertex->set_z(bbcz);
    vertex->set_z_err(0.6);
    vertex->set_t_err(m_tres);

    for (int iarm = 0; iarm < 2; iarm++)
    {
      vertex->set_bbc_ns(iarm, bbcn[iarm], bbcq[iarm], bbct[iarm]);
    }

    if (Verbosity() > 0)
    {
      std::cout << "bbc vertex z and t0 " << bbcz << ", "
                << bbct0 << std::endl;
    }

    m_bbcvertexmap->insert(vertex.release());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int BbcReconstruction::End(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int BbcReconstruction::createNodes(PHCompositeNode *topNode)
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

  m_bbcvertexmap = findNode::getClass<BbcVertexMap>(bbcNode, "BbcVertexMap");
  if (!m_bbcvertexmap)
  {
    m_bbcvertexmap = new BbcVertexMapv1();
    PHIODataNode<PHObject> *VertexMapNode = new PHIODataNode<PHObject>(m_bbcvertexmap, "BbcVertexMap", "PHObject");
    bbcNode->addNode(VertexMapNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
int BbcReconstruction::getNodes(PHCompositeNode *topNode)
{
  // BbcPmtContainer
  m_bbcpmts = findNode::getClass<BbcPmtContainer>(topNode, "BbcPmtContainer");
  if (!m_bbcpmts)
  {
    std::cout << PHWHERE << " BbcPmtContainer node not found on node tree" << std::endl;
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
