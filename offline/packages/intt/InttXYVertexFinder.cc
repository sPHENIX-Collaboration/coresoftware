#include "InttXYVertexFinder.h"
#include "INTTXYvtx.h"
#include "InttVertexMapv1.h"
#include "InttVertexv1.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/InttDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <cmath>
#include <string>

std::string xyinit_sDataType = "data";  // or "MC"
std::string xyinit_input_directory = "./";
std::string xyinit_file_name = "MC_ZF_xyvtx";
std::string xyinit_out_folder_directory = xyinit_input_directory + "/SemiFinal_RunXY_" + xyinit_file_name;

// std::pair<double, double> xyinit_beam_origin = {-0.015, 0.012};

std::pair<double, double> xyinit_beam_origin = {0.0, 0.0};
std::pair<double, double> xyinit_DCA_cut = {-1, 1};  // note : used
double xyinit_phi_diff_cut = 0.35;                   // note : used
int xyinit_clu_sum_adc_cut = -1;                     // note : used
int xyinit_clu_size_cut = 4;                         // note : used
int xyinit_N_clu_cut = 350;                          // note : used
int xyinit_N_clu_cutl = 20;                          // note : used

int xyinit_zvtx_cal_require = 15;
double xyinit_peek = 3.324;

bool GetPhiCheckTag(
    std::vector<INTTXYvtx::clu_info>& temp_sPH_inner_nocolumn_vec,
    std::vector<INTTXYvtx::clu_info>& temp_sPH_outer_nocolumn_vec)
{
  int inner_1_check = 0;
  int inner_2_check = 0;
  int inner_3_check = 0;
  int inner_4_check = 0;
  // for (unsigned int inner_i = 0; inner_i < temp_sPH_inner_nocolumn_vec.size(); inner_i++)
  for (const auto& clu : temp_sPH_inner_nocolumn_vec)
  {
    if (clu.phi >= 0 && clu.phi < 90)
    {
      inner_1_check = 1;
    }
    if (clu.phi >= 90 && clu.phi < 180)
    {
      inner_2_check = 1;
    }
    if (clu.phi >= 180 && clu.phi < 270)
    {
      inner_3_check = 1;
    }
    if (clu.phi >= 270 && clu.phi < 360)
    {
      inner_4_check = 1;
    }

    if ((inner_1_check + inner_2_check + inner_3_check + inner_4_check) == 4)
    {
      break;
    }
  }

  int outer_1_check = 0;
  int outer_2_check = 0;
  int outer_3_check = 0;
  int outer_4_check = 0;
  // for (unsigned int outer_i = 0; outer_i < temp_sPH_outer_nocolumn_vec.size(); outer_i++)
  for (const auto& clu : temp_sPH_outer_nocolumn_vec)
  {
    if (clu.phi >= 0 && clu.phi < 90)
    {
      outer_1_check = 1;
    }
    if (clu.phi >= 90 && clu.phi < 180)
    {
      outer_2_check = 1;
    }
    if (clu.phi >= 180 && clu.phi < 270)
    {
      outer_3_check = 1;
    }
    if (clu.phi >= 270 && clu.phi < 360)
    {
      outer_4_check = 1;
    }

    if ((outer_1_check + outer_2_check + outer_3_check + outer_4_check) == 4)
    {
      break;
    }
  }

  if ((inner_1_check + inner_2_check + inner_3_check + inner_4_check +
       outer_1_check + outer_2_check + outer_3_check + outer_4_check) != 8)
  {
    return false;
  }

  return true;
}

//____________________________________________________________________________..
InttXYVertexFinder::InttXYVertexFinder(const std::string& name)
  : SubsysReco(name)
  , m_inttxyvtx(new INTTXYvtx(
        xyinit_sDataType, xyinit_out_folder_directory, xyinit_beam_origin, xyinit_phi_diff_cut, xyinit_DCA_cut, xyinit_N_clu_cutl, xyinit_N_clu_cut, xyinit_peek))
{
  std::cout << "InttXYVertexFinder::InttXYVertexFinder(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
InttXYVertexFinder::~InttXYVertexFinder()
{
  std::cout << "InttXYVertexFinder::~InttXYVertexFinder() Calling dtor" << std::endl;
  delete m_inttxyvtx;
}

//____________________________________________________________________________..
int InttXYVertexFinder::Init(PHCompositeNode* /*topNode*/)
{
  std::cout << "InttXYVertexFinder::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InttXYVertexFinder::InitRun(PHCompositeNode* topNode)
{
  std::cout << "InttXYVertexFinder::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;

  // should be reset before init
  m_inttxyvtx->PrintMessageOpt(Verbosity() > 5);
  m_inttxyvtx->Init();

  if (createNodes(topNode) == Fun4AllReturnCodes::ABORTEVENT)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InttXYVertexFinder::process_event(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "InttXYVertexFinder::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  }

  static int event_i = 0;

  ActsGeometry* tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tGeometry)
  {
    std::cout << PHWHERE << "No ActsGeometry on node tree. Bailing." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  TrkrClusterContainer* clusterMap =
      findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!clusterMap)
  {
    std::cout << PHWHERE << "TrkrClusterContainer node is missing." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /////////////////////////
  // cluster list
  //  fill cluster
  std::vector<INTTXYvtx::clu_info> temp_sPH_inner_nocolumn_vec;
  std::vector<INTTXYvtx::clu_info> temp_sPH_outer_nocolumn_vec;

  std::vector<std::vector<double>> temp_sPH_nocolumn_vec(2);
  std::vector<std::vector<double>> temp_sPH_nocolumn_rz_vec(2);

  for (unsigned int inttlayer = 0; inttlayer < 4; inttlayer++)
  {
    std::vector<INTTXYvtx::clu_info>* p_temp_sPH_nocolumn_vec =
        (inttlayer < 2) ? (&temp_sPH_inner_nocolumn_vec)
                        : (&temp_sPH_outer_nocolumn_vec);

    for (const auto& hitsetkey : clusterMap->getHitSetKeys(TrkrDefs::TrkrId::inttId, inttlayer + 3))
    {
      auto range = clusterMap->getClusters(hitsetkey);

      for (auto clusIter = range.first; clusIter != range.second; ++clusIter)
      {
        const auto cluskey = clusIter->first;
        auto* const cluster = clusIter->second;

        const auto globalPos = tGeometry->getGlobalPosition(cluskey, cluster);
        double clu_x = globalPos.x() * 10.;  // convert to "mm" unit
        double clu_y = globalPos.y() * 10.;  // convert to "mm" unit
        double clu_z = globalPos.z() * 10.;  // convert to "mm" unit

        double clu_phi = (clu_y < 0) ? atan2(clu_y, clu_x) * (180. / M_PI) + 360
                                     : atan2(clu_y, clu_x) * (180. / M_PI);

        double clu_radius = sqrt(pow(clu_x, 2) + pow(clu_y, 2));

        p_temp_sPH_nocolumn_vec->push_back({-1,
                                            -1,
                                            (int) cluster->getAdc(),
                                            (int) cluster->getAdc(),
                                            (int) cluster->getPhiSize(),
                                            clu_x,                // change to "mm"
                                            clu_y,                // change to "mm"
                                            clu_z,                // change to "mm"
                                            (int) inttlayer + 3,  // note : should be 3 or 4, for the inner layer, this is for the mega cluster search
                                            clu_phi});

        temp_sPH_nocolumn_vec[0].push_back(clu_x);
        temp_sPH_nocolumn_vec[1].push_back(clu_y);

        temp_sPH_nocolumn_rz_vec[0].push_back(clu_z);
        temp_sPH_nocolumn_rz_vec[1].push_back((clu_phi > 180) ? clu_radius * -1 : clu_radius);
      }
    }
  }

  ////////////////////////
  //
  int NvtxMC = 1;
  double TrigZvtxMC = 0.;
  Long64_t bcofull = -1;

  m_inttxyvtx->ProcessEvt(
      event_i,
      temp_sPH_inner_nocolumn_vec,
      temp_sPH_outer_nocolumn_vec,
      temp_sPH_nocolumn_vec,
      temp_sPH_nocolumn_rz_vec,
      NvtxMC,
      TrigZvtxMC,
      GetPhiCheckTag(temp_sPH_inner_nocolumn_vec, temp_sPH_outer_nocolumn_vec),
      bcofull  // note : no bco_full for MC
  );

  m_inttxyvtx->ClearEvt();

  ///////////////////////////////////////////////////
  // calculate XY vertex
  if ((event_i % m_period) == 0)
  {
    // quadorant method
    std::vector<std::pair<double, double>> out_vtx = m_inttxyvtx->MacroVTXSquare(4, 10);
    m_vertex_quad[0] = out_vtx[0].first;
    m_vertex_quad[1] = out_vtx[0].second;

    if (Verbosity() > 1)
    {
      std::cout << " " << std::endl;
      std::cout << "The best vertex throughout the scan: " << out_vtx[0].first << " " << out_vtx[0].second << std::endl;
      std::cout << "The origin during that scan: " << out_vtx[1].first << " " << out_vtx[1].second << std::endl;
      std::cout << "Fit error, DCA and angle diff: " << out_vtx[2].first << " " << out_vtx[2].second << std::endl;
      std::cout << "fit pol0 pos Y, DCA and angle diff: " << out_vtx[3].first << " " << out_vtx[3].second << std::endl;
    }

    // line filled method
    std::vector<std::pair<double, double>> out_vtx_line = m_inttxyvtx->FillLine_FindVertex(
        {(out_vtx[0].first + out_vtx[1].first) / 2.,
         (out_vtx[0].second + out_vtx[1].second) / 2.},
        0.001);

    m_vertex_line[0] = out_vtx_line[0].first;
    m_vertex_line[1] = out_vtx_line[0].second;

    if (Verbosity() > 1)
    {
      std::cout << " " << std::endl;
      std::cout << "By fill-line method," << std::endl;
      std::cout << "Reco Run Vertex XY: " << out_vtx_line[0].first << " " << out_vtx_line[0].second << std::endl;
      std::cout << "Reco Run Vertex XY Error: " << out_vtx_line[1].first << " " << out_vtx_line[1].second << std::endl;
    }

    ///////////////////////////////////////////////////
    // convert xy-vertex from "mm" to "cm"  unit
    m_vertex_quad[0] *= 0.1;
    m_vertex_quad[1] *= 0.1;
    m_vertex_line[0] *= 0.1;
    m_vertex_line[1] *= 0.1;
  }

  auto vertex_quad = std::make_unique<InttVertexv1>();
  vertex_quad->set_x(m_vertex_quad[0]);
  vertex_quad->set_y(m_vertex_quad[1]);
  vertex_quad->set_z(m_vertex_quad[2]);
  // vertex_quad->set_id(0); // automatically set in insert
  //  vertex->set_z_err(); // no value asigned yet
  m_inttvertexmap->insert(vertex_quad.release());

  auto vertex_line = std::make_unique<InttVertexv1>();
  vertex_line->set_x(m_vertex_line[0]);
  vertex_line->set_y(m_vertex_line[1]);
  vertex_line->set_z(m_vertex_line[2]);
  // vertex_line->set_id(1); // automatically set in insert
  //  vertex->set_z_err(); // no value asigned yet
  m_inttvertexmap->insert(vertex_line.release());

  if (
      Verbosity() > 0 &&
      (event_i % m_period) == 0)
  {
    std::cout << event_i << ", vtx_xy_quad: " << m_vertex_quad[0] << " " << m_vertex_quad[1] << " " << m_vertex_quad[2]
              << ", vtx_xy_line: " << m_vertex_line[0] << " " << m_vertex_line[1] << " " << m_vertex_line[2]
              << std::endl;
  }
  event_i++;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InttXYVertexFinder::End(PHCompositeNode* /*topNode*/)
{
  m_inttxyvtx->PrintPlots();
  m_inttxyvtx->EndRun();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void InttXYVertexFinder::Print(const std::string& what) const
{
  std::cout << "InttXYVertexFinder::Print(const std::string &what) const Printing info for " << what << std::endl;
}

int InttXYVertexFinder::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // ---
  PHCompositeNode* inttNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "INTT"));
  if (!inttNode)
  {
    inttNode = new PHCompositeNode("INTT");
    dstNode->addNode(inttNode);
  }

  m_inttvertexmap = findNode::getClass<InttVertexMap>(inttNode, "InttVertexMap");
  if (!m_inttvertexmap)
  {
    m_inttvertexmap = new InttVertexMapv1();
    PHIODataNode<PHObject>* VertexMapNode = new PHIODataNode<PHObject>(m_inttvertexmap, "InttVertexMap", "PHObject");
    inttNode->addNode(VertexMapNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void InttXYVertexFinder::SetBeamCenter(const double beamx, const double beamy)
{
  if (m_inttxyvtx != nullptr)
  {
    m_inttxyvtx->SetBeamOrigin(beamx * 10., beamy * 10.);  // converted from "cm" to "mm"
  }
}

void InttXYVertexFinder::SetSaveHisto(const bool savehist)
{
  if (m_inttxyvtx != nullptr)
  {
    m_inttxyvtx->SetSaveHisto(savehist);
  }
}

void InttXYVertexFinder::EnableDrawHisto(const bool enable)
{
  if (m_inttxyvtx != nullptr)
  {
    m_inttxyvtx->EnableDrawHisto(enable);
  }
}

void InttXYVertexFinder::EnableQA(const bool enable)
{
  if (m_inttxyvtx != nullptr)
  {
    m_inttxyvtx->EnableQA(enable);
  }
}
