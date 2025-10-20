#include "InttZVertexFinder.h"
#include "INTTZvtx.h"
#include "InttVertexMapv1.h"
#include "InttVertexv1.h"

/// Fun4All includes
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
#include <filesystem>

/////////////////////////////////////////
// init value for temporal use
std::string sDataType = "data";  // or "MC"
std::string sOutFolderDirectory = "./";
std::pair<double, double> beam_origin = {-0.23436750, 2.5985125};  // note : for run20869
double phi_diff_cut = 1;                                           // 0.265 + 0.269;            // note : for run20869
// std::pair<double, double> DCA_cut             = {0.277 - 0.730, 0.277 + 0.730}; // note : for run20869
std::pair<double, double> DCA_cut = {-3, 3};  // note : for run20869
int N_clu_cut = 10000;
int N_clu_cutl = 4;                 // 20
unsigned int zvtx_cal_require = 3;  // 5; //15;
std::pair<double, double> zvtx_QA_width = {40, 70};
double zvtx_QA_ratio = 0.00006;
bool draw_event_display = false;
bool enable_qa = false;

int clu_sum_adc_cut = 31;
int clu_size_cut = 4;

int data_type = 2;  // note : private gen data cluster

/*
bool GetPhiCheckTag(
        std::vector<INTTZvtx::clu_info>& temp_sPH_inner_nocolumn_vec,
        std::vector<INTTZvtx::clu_info>& temp_sPH_outer_nocolumn_vec
      )
{
    int inner_1_check = 0;
    int inner_2_check = 0;
    int inner_3_check = 0;
    int inner_4_check = 0;
    for (unsigned int inner_i = 0; inner_i < temp_sPH_inner_nocolumn_vec.size(); inner_i++) {
        if (temp_sPH_inner_nocolumn_vec[inner_i].phi >=   0 && temp_sPH_inner_nocolumn_vec[inner_i].phi <  90) inner_1_check = 1;
        if (temp_sPH_inner_nocolumn_vec[inner_i].phi >=  90 && temp_sPH_inner_nocolumn_vec[inner_i].phi < 180) inner_2_check = 1;
        if (temp_sPH_inner_nocolumn_vec[inner_i].phi >= 180 && temp_sPH_inner_nocolumn_vec[inner_i].phi < 270) inner_3_check = 1;
        if (temp_sPH_inner_nocolumn_vec[inner_i].phi >= 270 && temp_sPH_inner_nocolumn_vec[inner_i].phi < 360) inner_4_check = 1;

        if ( (inner_1_check + inner_2_check + inner_3_check + inner_4_check) == 4 ) break;
    }

    int outer_1_check = 0;
    int outer_2_check = 0;
    int outer_3_check = 0;
    int outer_4_check = 0;
    for (unsigned int outer_i = 0; outer_i < temp_sPH_outer_nocolumn_vec.size(); outer_i++) {
        if (temp_sPH_outer_nocolumn_vec[outer_i].phi >=   0 && temp_sPH_outer_nocolumn_vec[outer_i].phi <  90) outer_1_check = 1;
        if (temp_sPH_outer_nocolumn_vec[outer_i].phi >=  90 && temp_sPH_outer_nocolumn_vec[outer_i].phi < 180) outer_2_check = 1;
        if (temp_sPH_outer_nocolumn_vec[outer_i].phi >= 180 && temp_sPH_outer_nocolumn_vec[outer_i].phi < 270) outer_3_check = 1;
        if (temp_sPH_outer_nocolumn_vec[outer_i].phi >= 270 && temp_sPH_outer_nocolumn_vec[outer_i].phi < 360) outer_4_check = 1;

        if ( (outer_1_check + outer_2_check + outer_3_check + outer_4_check) == 4 ) break;
    }

    if ( (inner_1_check + inner_2_check + inner_3_check + inner_4_check +
          outer_1_check + outer_2_check + outer_3_check + outer_4_check) != 8 ) {return false;}
    else { return true; }
}
*/

//____________________________________________________________________________..
InttZVertexFinder::InttZVertexFinder(const std::string& name)
  : SubsysReco(name)
  , m_inttzvtx(new INTTZvtx(sDataType,
                            sOutFolderDirectory,
                            beam_origin,
                            phi_diff_cut,
                            DCA_cut,
                            N_clu_cutl,
                            N_clu_cut,
                            zvtx_cal_require,
                            zvtx_QA_width,
                            draw_event_display,
                            enable_qa))
{
  std::cout << "InttZVertexFinder::InttZVertexFinder(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
InttZVertexFinder::~InttZVertexFinder()
{
  std::cout << "InttZVertexFinder::~InttZVertexFinder() Calling dtor" << std::endl;
  delete m_inttzvtx;
}

//____________________________________________________________________________..
int InttZVertexFinder::Init(PHCompositeNode* /*topNode*/)
{
  std::cout << "InttZVertexFinder::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  // m_inttzvtx->Init();
  //
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InttZVertexFinder::InitRun(PHCompositeNode* topNode)
{
  std::cout << "InttZVertexFinder::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;

  if (createNodes(topNode) == Fun4AllReturnCodes::ABORTEVENT)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_inttzvtx->SetPrintMessageOpt(Verbosity() > 0);

  m_inttzvtx->Init();

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttZVertexFinder::createNodes(PHCompositeNode* topNode)
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

//____________________________________________________________________________..
int InttZVertexFinder::process_event(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "InttZVertexFinder::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  }

  //  INTTClu -> EvtInit(event_i);
  //  INTTClu -> EvtSetCluGroup();

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
  std::vector<INTTZvtx::clu_info> temp_sPH_inner_nocolumn_vec;
  std::vector<INTTZvtx::clu_info> temp_sPH_outer_nocolumn_vec;

  std::vector<std::vector<double> > temp_sPH_nocolumn_vec(2);
  std::vector<std::vector<double> > temp_sPH_nocolumn_rz_vec(2);

  for (unsigned int inttlayer = 0; inttlayer < 4; inttlayer++)
  {
    std::vector<INTTZvtx::clu_info>* p_temp_sPH_nocolumn_vec =
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
        double clu_x = globalPos.x() * 10.;  // convert from "cm" to "mm" unit
        double clu_y = globalPos.y() * 10.;  // convert from "cm" to "mm" unit
        double clu_z = globalPos.z() * 10.;  // convert from "cm" to "mm" unit

        double clu_phi = (clu_y < 0) ? atan2(clu_y, clu_x) * (180. / M_PI) + 360
                                     : atan2(clu_y, clu_x) * (180. / M_PI);

        double clu_radius = sqrt(pow(clu_x, 2) + pow(clu_y, 2));

        // int size = cluster->getSize();
        int adc = cluster->getAdc();

        if (adc < 40)
        {
          continue;
        }

        p_temp_sPH_nocolumn_vec->push_back({-1,
                                            -1,
                                            (int) cluster->getAdc(),
                                            (int) cluster->getAdc(),
                                            (int) cluster->getPhiSize(),
                                            clu_x,
                                            clu_y,
                                            clu_z,
                                            (int) inttlayer + 3,  // note : should be 3 or 4, for the inner layer, this is for the mega cluster search
                                            clu_phi});

        temp_sPH_nocolumn_vec[0].push_back(clu_x);
        temp_sPH_nocolumn_vec[1].push_back(clu_y);

        temp_sPH_nocolumn_rz_vec[0].push_back(clu_z);
        temp_sPH_nocolumn_rz_vec[1].push_back((clu_phi > 180) ? clu_radius * -1 : clu_radius);

        // std::cout<<"cluster pos: ";
        // std::cout<<globalPos.x()<<" ";
        // std::cout<<globalPos.y()<<" ";
        // std::cout<<globalPos.z()<<std::endl;
      }
    }
  }

  ////////////////////////

  int NvtxMC = 1;
  double TrigZvtxMC = 0.;
  Long64_t bcofull = 0;

  bool status = m_inttzvtx->ProcessEvt(
      event_i,
      temp_sPH_inner_nocolumn_vec,
      temp_sPH_outer_nocolumn_vec,
      temp_sPH_nocolumn_vec,
      temp_sPH_nocolumn_rz_vec,
      NvtxMC,
      TrigZvtxMC,
      true,  // GetPhiCheckTag(temp_sPH_inner_nocolumn_vec, temp_sPH_outer_nocolumn_vec),
      bcofull,
      5  // centrality bin, note : no bco_full for MC
  );

  std::cout << "InttZVertex:process_evt status = " << (status ? "good" : "failed") << std::endl;

  std::vector<double> vtxout = m_inttzvtx->GetEvtZPeak();              // mm unit
  std::pair<double, double> beamorigin = m_inttzvtx->GetBeamOrigin();  // mm unit

  std::vector<double> final_vtx{beamorigin.first * 0.1,
                                beamorigin.second * 0.1,
                                vtxout[1] * 0.1};  // mm -> cm by 0.1

  INTTZvtx::ZvtxInfo& zvtxinfo = m_inttzvtx->GetZvtxInfo();

  //--if(m_inttvertexmap->size()>0){
  //--  std::cout<<"added to existing vertex"<<std::endl;
  //--  for(auto vtx = m_inttvertexmap->begin(); vtx!=m_inttvertexmap->end(); ++vtx){
  //--    vtx->second->set_x(final_vtx[0]);
  //--    vtx->second->set_y(final_vtx[1]);
  //--    vtx->second->set_z(final_vtx[2]);
  //--  }
  //--}
  //--else {
  //--std::cout<<"new vertex created"<<std::endl;
  auto vertex = std::make_unique<InttVertexv1>();
  vertex->set_x(final_vtx[0]);
  vertex->set_y(final_vtx[1]);
  vertex->set_z(final_vtx[2]);

  vertex->set_chi2ndf(zvtxinfo.chi2ndf);
  vertex->set_width(zvtxinfo.width);
  vertex->set_good(zvtxinfo.good);
  vertex->set_nclus(zvtxinfo.nclus);
  vertex->set_ntracklet(zvtxinfo.ntracklets);
  vertex->set_ngroup(zvtxinfo.ngroup);
  vertex->set_peakratio(zvtxinfo.peakratio);
  vertex->set_peakwidth(zvtxinfo.peakwidth);

  // vertex->set_z_err(); // no value asigned yet
  m_inttvertexmap->insert(vertex.release());
  //}

  if (Verbosity() > 0)
  {
    std::cout << "vecsize in : " << temp_sPH_inner_nocolumn_vec.size()
              << ", out: " << temp_sPH_outer_nocolumn_vec.size()
              << ", zvtx : " << vtxout[1] << " +- " << vtxout[2]
              << " sts : " << vtxout[0]
              << std::endl;
  }

  m_inttzvtx->ClearEvt();

  event_i++;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InttZVertexFinder::End(PHCompositeNode* /*topNode*/)
{
  if (Verbosity())
  {
    std::cout << "InttZVertexFinder::End(PHCompositeNode *topNode) " << std::endl;
  }

  m_inttzvtx->PrintPlots();
  m_inttzvtx->EndRun();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void InttZVertexFinder::Print(const std::string& what) const
{
  std::cout << "InttZVertexFinder::Print(const std::string &what) const Printing info for " << what << std::endl;
}

void InttZVertexFinder::SetBeamCenter(const double beamx, const double beamy)
{
  if (m_inttzvtx != nullptr)
  {
    m_inttzvtx->SetBeamOrigin(beamx * 10., beamy * 10.);  // convert to cm to mm unit
  }
}

void InttZVertexFinder::SetOutDirectory(const std::string& outDirectory)
{
  if (m_inttzvtx != nullptr)
  {
    m_inttzvtx->SetOutDirectory(outDirectory);  // convert to cm to mm unit
  }
}

void InttZVertexFinder::EnableQA(const bool enableQA)
{
  if (m_inttzvtx != nullptr)
  {
    m_inttzvtx->EnableQA(enableQA);
  }
}

void InttZVertexFinder::EnableEventDisplay(const bool enableEvtDisp)
{
  if (m_inttzvtx != nullptr)
  {
    m_inttzvtx->EnableEventDisplay(enableEvtDisp);
  }
}
