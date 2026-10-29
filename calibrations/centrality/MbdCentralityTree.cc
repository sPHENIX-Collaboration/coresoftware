#include "MbdCentralityTree.h"

#include <ffarawobjects/Gl1Packet.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <zdcinfo/Zdcinfo.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <iostream>

MbdCentralityTree::MbdCentralityTree(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , m_filename(filename)
{
}

int MbdCentralityTree::Init(PHCompositeNode * /*topNode*/)
{
  m_file = new TFile(m_filename.c_str(), "RECREATE");
  m_tree = new TTree("T", "MBD and ZDC per event, for the centrality calibration");
  m_tree->Branch("gl1_scaled", &m_gl1_scaled, "gl1_scaled/l");
  m_tree->Branch("mbd_vertex", &m_mbd_vertex, "mbd_vertex/F");
  m_tree->Branch("mbd_time_zero", &m_mbd_time_zero, "mbd_time_zero/F");
  m_tree->Branch("mbd_charge", m_mbd_charge, "mbd_charge[128]/F");
  m_tree->Branch("mbd_time", m_mbd_time, "mbd_time[128]/F");
  m_tree->Branch("zdc_sum", m_zdc_sum, "zdc_sum[2]/F");
  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdCentralityTree::process_event(PHCompositeNode *topNode)
{
  Gl1Packet *gl1 = findNode::getClass<Gl1Packet>(topNode, 14001);
  MbdPmtContainer *pmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  Zdcinfo *zdcinfo = findNode::getClass<Zdcinfo>(topNode, "Zdcinfo");
  if (!gl1 || !pmts || !vertexmap || !zdcinfo)
  {
    std::cout << PHWHERE << " missing node:" << (gl1 ? "" : " 14001 (Gl1Packet)") << (pmts ? "" : " MbdPmtContainer")
              << (vertexmap ? "" : " GlobalVertexMap") << (zdcinfo ? "" : " Zdcinfo") << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_gl1_scaled = gl1->getScaledVector();
  bool triggered = false;
  for (int bit : m_trigger_bits)
  {
    triggered = triggered || ((m_gl1_scaled >> bit) & 0x1U);
  }
  if (!triggered)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // 999: no MBD vertex in this event; -999: a vertex without a valid z
  m_mbd_vertex = 999;
  m_mbd_time_zero = 999;
  std::vector<GlobalVertex *> vertices = vertexmap->get_gvtxs_with_type({GlobalVertex::VTXTYPE::MBD});
  if (!vertices.empty() && vertices.at(0))
  {
    m_mbd_vertex = vertices.at(0)->get_z();
    m_mbd_time_zero = vertices.at(0)->get_t();
  }
  if (std::isnan(m_mbd_vertex))
  {
    m_mbd_vertex = -999;
  }

  for (int ipmt = 0; ipmt < 128; ipmt++)
  {
    MbdPmtHit *hit = pmts->get_pmt(ipmt);
    m_mbd_charge[ipmt] = hit->get_q();
    m_mbd_time[ipmt] = hit->get_time();
  }
  for (int iside = 0; iside < 2; iside++)
  {
    m_zdc_sum[iside] = zdcinfo->get_zdc_energy(iside);
  }

  m_tree->Fill();
  m_nevents++;
  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdCentralityTree::End(PHCompositeNode * /*topNode*/)
{
  std::cout << Name() << ": " << m_nevents << " events written to " << m_filename << std::endl;
  m_file->cd();
  m_tree->Write();
  m_file->Close();
  delete m_file;
  m_file = nullptr;
  return Fun4AllReturnCodes::EVENT_OK;
}
