#include "InttVertexFinder.h"
#include "InttVertexMapv1.h"
#include "InttVertexv1.h"

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase/ActsGeometry.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                         // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TH1.h>

#include <cmath>
#include <iostream>
#include <map>       // for _Rb_tree_const_i...
#include <memory>    // for make_unique, uni...
#include <utility>   // for move, pair
#include <vector>    // for vector

InttVertexFinder::InttVertexFinder(const std::string& name)
  : SubsysReco(name)
{}

InttVertexFinder::~InttVertexFinder()
{
  delete h_zvtxseed_;
}

int InttVertexFinder::Init(PHCompositeNode* /*topNode*/)
{

  delete h_zvtxseed_;
  h_zvtxseed_ = new TH1F("h_zvtxseed", "Zvertex Seed histogram", 200, -50, 50);

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttVertexFinder::InitRun(PHCompositeNode* topNode)
{

  if( createNodes(topNode) == Fun4AllReturnCodes::ABORTEVENT) {
     return Fun4AllReturnCodes::ABORTEVENT;
  }

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    std::cout << "====================== InttVertexFinder::InitRun() =====================" << std::endl;
    std::cout<<"            beamcenter "<<xbeam_<<" "<<ybeam_<<std::endl;
    std::cout << "===========================================================================" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttVertexFinder::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // ---
  PHCompositeNode *inttNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "INTT"));
  if (!inttNode)
  {
    inttNode = new PHCompositeNode("INTT");
    dstNode->addNode(inttNode);
  }

  m_inttvertexmap = findNode::getClass<InttVertexMap>(inttNode, "InttVertexMap");
  if (!m_inttvertexmap)
  {
    m_inttvertexmap = new InttVertexMapv1();
    PHIODataNode<PHObject> *VertexMapNode = new PHIODataNode<PHObject>(m_inttvertexmap, "InttVertexMap", "PHObject");
    inttNode->addNode(VertexMapNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}



int InttVertexFinder::process_event(PHCompositeNode* topNode)
{

  if (Verbosity() > 5)
  {
    std::cout << "Beginning process_event in InttVertexFinder" << std::endl;
  }
  
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << "No ActsGeometry on node tree. Bailing." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  
  
  // get node for clusters
  m_clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterlist)
  {
    std::cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTER." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

/////////////////////////////////////////

  double zcenter, zrms, zmean;
  double zvertex = calculateZvertex(&zcenter, &zrms, &zmean);


  auto vertex = std::make_unique<InttVertexv1>();
  vertex->set_z(zvertex);
  //vertex->set_z_err(); // no value asigned yet

  if (Verbosity() > 0)
  {
    std::cout << "intt vertex z " << zvertex << std::endl;
  }

  m_inttvertexmap->insert(vertex.release());


///////////////////
  return Fun4AllReturnCodes::EVENT_OK;
}

double InttVertexFinder::calculateZvertex(
       double*  zcenter, double* zrms, double* zmean
       )
{
  struct ClustInfo {
    int layer;
    int adc;
    Acts::Vector3 pos;
  };


  int nCluster=0;
  //bool exceedNwrite=false;
  std::vector<ClustInfo> clusters[8][2]; // phi_angle:0-7, layer: inner=0; outer=1

  for (unsigned int inttlayer = 0; inttlayer < 4; inttlayer++)
  {
    int inout = (inttlayer<2) ? 0 : 1;
    for (const auto &hitsetkey : m_clusterlist->getHitSetKeys(TrkrDefs::TrkrId::inttId, inttlayer + 3))
    {
      auto range = m_clusterlist->getClusters(hitsetkey);

      for (auto clusIter = range.first; clusIter != range.second; ++clusIter)
      {
        const auto cluskey = clusIter->first;
        const auto cluster = clusIter->second;
        const auto globalPos = m_tGeometry->getGlobalPosition(cluskey, cluster);
        //int ladder_z   = InttDefs::getLadderZId(cluskey);
        //int ladder_phi = InttDefs::getLadderPhiId(cluskey);
        //int size       = cluster->getSize();

        //--if(nCluster<5) {
        //--  std::cout<<"xyz : "<<globalPos.x()<<" "<< globalPos.y()<<" "<< globalPos.z()<<" :  "
        //--    <<cluster->getAdc()<<" "<<size<<" "<<inttlayer<<" "<<ladder_z<<" "<<ladder_phi<<std::endl;
        //--}
        //--else {
        //--  if(!exceedNwrite) {
        //--    std::cout<<" exceed : ncluster limit.  no more cluster xyz printed"<<std::endl;
        //--    exceedNwrite=true;
        //--  }
        //--}

        ClustInfo info;
        info.layer = inttlayer;
        info.adc   = cluster->getAdc();
        info.pos   = globalPos;

	double phi = atan2(globalPos.y(), globalPos.x());

	int iphi = (phi + M_PI)/M_PI/4.;
	if(iphi<0)  iphi+=8;
	iphi%=8;
	//std::cout<<"phi : "<<phi<<" "<<iphi<<std::endl;


	clusters[iphi][inout].push_back(info);

        nCluster++;
      }
    }
  }


  Acts::Vector3 beamspot(xbeam_, ybeam_, 0);
  std::vector<double> vz_array;

  for(int iphi=0; iphi<8; iphi++){
    for(auto c1=clusters[iphi][0].begin(); c1!=clusters[iphi][0].end(); ++c1) // inner
    {
      for(auto c2=clusters[iphi][1].begin(); c2!=clusters[iphi][1].end(); ++c2) // outer
      {
        if(c1->adc<40 || c2->adc<40) continue;

        //TVector3 p1 = c1->pos - beamspot;
        //TVector3 p2 = c2->pos - beamspot;
        Acts::Vector3 p1 = c1->pos - beamspot;
        Acts::Vector3 p2 = c2->pos - beamspot;
        // skip bad compbination
        double p1_ang = atan2(p1.y(), p1.x());
        double p2_ang = atan2(p2.y(), p2.x());
        double d_ang = p2_ang - p1_ang;

        if(fabs(d_ang)>0.2 ) continue;

        //TVector3 u = p2 - p1;
        Acts::Vector3 u = p2 - p1;
        double unorm = sqrt(u.x()*u.x()+u.y()*u.y());
        if(unorm<0.00001) continue;



        // unit vector in 2D
        double ux = u.x()/unorm;
        double uy = u.y()/unorm;
        double uz = u.z()/unorm; // same normalization factor(2D) is multiplied

        Acts::Vector3 p0   = beamspot - p1;
        //TVector3 p0   = beamspot - p1;

        double dca_p0 = p0.x()*uy - p0.y()*ux; // cross product of p0 x u
        double len_p0 = p0.x()*ux + p0.y()*uy; // dot product of p0 . u


        // beam center in X-Y plane
        //double vx = len_p0*ux + p1.x();
        //double vy = len_p0*uy + p1.y();

        double vz = len_p0*uz + p1.z();

        if(fabs(d_ang)<0.05 && fabs(dca_p0) < 1.0)
        {
          h_zvtxseed_->Fill(vz);
          vz_array.push_back(vz);
        }

      }
    }
  }

  // calculate trancated mean of DCA~Z histogram as Z-vertex position
  double zvtx=-9999.;

  double zcenter1 = -9999.;
  double zmean1   = -9999.;
  double zrms1    = -9999.;
  if(vz_array.size()>3){
    double zbin    = h_zvtxseed_->GetMaximumBin();
    zcenter1       = h_zvtxseed_->GetBinCenter(zbin);
    zmean1         = h_zvtxseed_->GetMean();
    zrms1          = h_zvtxseed_->GetRMS();
    if(zrms1<20) zrms1 = 20;

    double zmax = zcenter1 + zrms1; // 1 sigma
    double zmin = zcenter1 - zrms1; // 1 sigma

    double zsum=0.;
    int    zcount=0;
    for(auto iz=vz_array.begin(); iz!=vz_array.end(); ++iz){
      double vz = (*iz);
      if(zmin<vz&&vz<zmax){
        zsum+=vz;
        zcount++;
      }
    }
    if(zcount>0) zvtx = zsum/zcount;

    if (Verbosity() > 0)
    {
      std::cout<<"ZVTX: "<<zvtx<<" "<<zcenter1<<" "<<zmean1<<" "<<zrms1<<" "<<zbin<<std::endl; //" "<<mbdt.bz<<std::endl;
    }
  }

  if(zcenter!=nullptr) *zcenter = zcenter1;
  if(zrms   !=nullptr) *zrms    = zrms1;
  if(zmean  !=nullptr) *zmean   = zmean1;

  return zvtx;
}


