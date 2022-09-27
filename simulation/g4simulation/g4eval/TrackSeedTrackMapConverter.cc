
#include "TrackSeedTrackMapConverter.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>

#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                        // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                      // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h> 

//____________________________________________________________________________..
TrackSeedTrackMapConverter::TrackSeedTrackMapConverter(const std::string &name):
 SubsysReco(name)
{
}

//____________________________________________________________________________..
TrackSeedTrackMapConverter::~TrackSeedTrackMapConverter()
{
}

//____________________________________________________________________________..
int TrackSeedTrackMapConverter::InitRun(PHCompositeNode *topNode)
{
  int ret = getNodes(topNode);
  return ret;
}

//____________________________________________________________________________..
int TrackSeedTrackMapConverter::process_event(PHCompositeNode*)
{
  if(Verbosity() > 1)
    {
      std::cout <<"silicon seed map size " << m_siContainer->size() << std::endl;
      for(auto iter = m_siContainer->begin(); iter != m_siContainer->end();
	  ++iter)
	{
	  auto seed = *iter;
	  if(!seed){
	    std::cout << "no seed at index "<< m_siContainer->index(iter) 
		      << std::endl;
	    continue;
	  }
	  seed->identify();
	}
      if(m_tpcContainer)
	{
	  std::cout << "tpc seed map size " << m_tpcContainer->size() << std::endl;
	  for(auto iter = m_tpcContainer->begin(); iter != m_tpcContainer->end();
	      ++iter)
	    {
	      auto seed = *iter;
	      if(!seed) {
		std::cout << "no tpc seed at entry " << m_tpcContainer->index(iter) 
			  << std::endl;
		continue;
	      }
	      seed->identify();
	    }
	}
    }

  unsigned int trackid = 0;
  for(const auto& trackSeed : *m_seedContainer)
    {
      /// If the seed was removed, skip it
      if(!trackSeed)
	{ continue; }

      if(m_trackSeedName.find("SvtxTrackSeed") != std::string::npos)
	{
	  /// Catches entries in the vector removed by ghost finder
	  unsigned int tpcindex = trackSeed->get_tpc_seed_index();
	  TrackSeed* seed = m_tpcContainer->get(tpcindex);
	  if(!seed)
	    { continue; }
	}

      auto svtxtrack = std::make_unique<SvtxTrack_v4>();

      if(Verbosity() > 0)
	{
	  std::cout << "iterating track seed " << trackid << std::endl;
	}

      svtxtrack->set_id(trackid);
      trackid++;

      /// If we've run the track matching
      if(m_trackSeedName.find("SvtxTrackSeed") != std::string::npos)
	{
	  if(Verbosity() > 0)
	    {
	      std::cout << "tpc seed id " << trackSeed->get_tpc_seed_index() <<std::endl;
	      std::cout << "si seed id " << trackSeed->get_silicon_seed_index() << std::endl;
	    }

	  unsigned int seedindex = trackSeed->get_tpc_seed_index();
	  TrackSeed *tpcseed = m_tpcContainer->get(seedindex);
	  if(trackSeed->get_silicon_seed_index() == std::numeric_limits<unsigned int>::max())
	    {      
	      /// Didn't find a match, so just use the tpc seed
	      svtxtrack->set_x(tpcseed->get_x());
	      svtxtrack->set_y(tpcseed->get_y());
	      svtxtrack->set_z(tpcseed->get_z()); 
	    }
	  else
	    {
	      TrackSeed *siseed = m_siContainer->get(trackSeed->get_silicon_seed_index());
	      svtxtrack->set_x(siseed->get_x());
	      svtxtrack->set_y(siseed->get_y());
	      svtxtrack->set_z(siseed->get_z());
	      addKeys(svtxtrack, siseed);
	      svtxtrack->set_silicon_seed(siseed);
	    }


	  svtxtrack->set_charge( tpcseed->get_qOverR() > 0 ? 1 : -1);
	  svtxtrack->set_px(tpcseed->get_px(m_clusters,m_tGeometry));
	  svtxtrack->set_py(tpcseed->get_py(m_clusters,m_tGeometry));
	  svtxtrack->set_pz(tpcseed->get_pz());
	  
	  addKeys(svtxtrack, tpcseed);
	  svtxtrack->set_tpc_seed(tpcseed);
	

	}
      else
	{
	  /// Otherwise we are using an individual subdetectors container

	  svtxtrack->set_x(trackSeed->get_x());
	  svtxtrack->set_y(trackSeed->get_y());
	  svtxtrack->set_z(trackSeed->get_z());
	  svtxtrack->set_charge( trackSeed->get_qOverR() > 0 ? 1 : -1);
	  svtxtrack->set_px(trackSeed->get_px(m_clusters,m_tGeometry));
	  svtxtrack->set_py(trackSeed->get_py(m_clusters,m_tGeometry));
	  svtxtrack->set_pz(trackSeed->get_pz());

          // calculate chisq and ndf
          double R = 1./fabs(trackSeed->get_qOverR());
          double X0 = trackSeed->get_X0();
          double Y0 = trackSeed->get_Y0();
          double Z0 = trackSeed->get_Z0();
          double slope = trackSeed->get_slope();
          std::vector<double> xy_error2;
          std::vector<double> rz_error2;
          std::vector<double> xy_residuals;
          std::vector<double> rz_residuals;
          std::vector<double> x_circle;
          std::vector<double> y_circle;
          std::vector<double> z_line;
          for(auto c_iter = trackSeed->begin_cluster_keys();
              c_iter != trackSeed->end_cluster_keys();
              ++c_iter)
          {
            TrkrCluster* c = m_clusters->findCluster(*c_iter);
            Acts::Vector3 pos = m_tGeometry->getGlobalPosition(*c_iter,c);
            double x = pos(0);
            double y = pos(1);
            double z = pos(2);
            double r = sqrt(x*x+y*y);
            double dx = x-X0;
            double dy = y-Y0;
            double xy_centerdist = sqrt(dx*dx+dy*dy);
            // method lifted from ALICEKF::GetCircleClusterResiduals
            xy_residuals.push_back(xy_centerdist-R);
            // method lifted from ALICEKF::GetLineClusterResiduals
            rz_residuals.push_back(fabs(-slope*r+z-Z0)/sqrt(slope*slope+1));

            // ignoring covariance for simplicity
            xy_error2.push_back(c->getActsLocalError(0,0)+c->getActsLocalError(1,1));
            rz_error2.push_back(c->getActsLocalError(2,2));
            double phi = atan2(dy,dx);
            x_circle.push_back(R*cos(phi)+X0);
            y_circle.push_back(R*sin(phi)+Y0);
            z_line.push_back(R*slope+Z0);
          }
          double chi2 = 0.;
          for(unsigned int i=0; i<xy_residuals.size(); i++)
          {
            if(std::isnan(xy_error2[i])) xy_error2[i] = 0.01;
            if(std::isnan(rz_error2[i])) rz_error2[i] = 0.01;
            // method lifted from GPUTPCTrackParam::Filter
            chi2 += xy_residuals[i]*xy_residuals[i]/xy_error2[i] + rz_residuals[i]*rz_residuals[i]/rz_error2[i];
          }
          svtxtrack->set_chisq(chi2);
          // GPUTPCTrackParam initially sets NDF to -3 on first cluster and increments by 2 with every application of filter
          svtxtrack->set_ndf(2*xy_residuals.size()-5);

          addKeys(svtxtrack, trackSeed);
	  if(m_trackSeedName.find("SiliconTrackSeed") != std::string::npos)
	    {
	      svtxtrack->set_silicon_seed(trackSeed);
	      svtxtrack->set_tpc_seed(nullptr);
	    }
	  else if(m_trackSeedName.find("TpcTrackSeed") != std::string::npos)
	    {
	      svtxtrack->set_tpc_seed(trackSeed);
	      svtxtrack->set_silicon_seed(nullptr);
	    }
	}

      if(Verbosity() > 0)
	{
	  std::cout << "Inserting svtxtrack into map " << std::endl;
	  svtxtrack->identify();
	}
      
      m_trackMap->insert(svtxtrack.get());

    }

  return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________..
int TrackSeedTrackMapConverter::End(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void TrackSeedTrackMapConverter::addKeys(std::unique_ptr<SvtxTrack_v4>& track, TrackSeed *seed)
{
  for(TrackSeed::ConstClusterKeyIter iter = seed->begin_cluster_keys();
      iter != seed->end_cluster_keys();
      ++iter)
    {
      track->insert_cluster_key(*iter);
    }
}


int TrackSeedTrackMapConverter::getNodes(PHCompositeNode *topNode)
{
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  if(!m_trackMap)
    {
      // create it
      PHNodeIterator iter(topNode);
      
      PHCompositeNode* dstNode = static_cast<PHCompositeNode*>(iter.findFirst(
									      "PHCompositeNode", "DST"));
      if (!dstNode)
	{
	  std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
      PHNodeIterator iter_dst(dstNode);
      
      // Create the SVTX node
      PHCompositeNode* tb_node =
	dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode",
							  "SVTX"));
      if (!tb_node)
	{
	  tb_node = new PHCompositeNode("SVTX");
	  dstNode->addNode(tb_node);
	  if (Verbosity() > 0)
	    { std::cout << PHWHERE << "SVTX node added" << std::endl; }
	}
      
      
      m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
      if (!m_trackMap)
	{
	  m_trackMap = new SvtxTrackMap_v1;
	  PHIODataNode<PHObject>* tracks_node = 
	    new PHIODataNode<PHObject>(m_trackMap, m_trackMapName, "PHObject");
	  tb_node->addNode(tracks_node);
	  if (Verbosity() > 0){
	    std::cout << PHWHERE << "Svtx/" << m_trackMapName  << " node added" << std::endl;
	  }
	}
    }
  
  m_seedContainer = findNode::getClass<TrackSeedContainer>(topNode, m_trackSeedName);
  if(!m_seedContainer)
    {
      std::cout << PHWHERE << " Can't find track seed container " << m_trackSeedName << ", can't continue."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_tpcContainer = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if(!m_tpcContainer)
    {
      std::cout << PHWHERE << "WARNING, TrackSeedTrackMapConverter may seg fault depending on what seeding algorithm this is run after" << std::endl;
    }

  m_siContainer = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if(!m_siContainer)
    {
        std::cout << PHWHERE << "WARNING, TrackSeedTrackMapConverter may seg fault depending on what seeding algorithm this is run after" << std::endl;
    }

  m_clusters = findNode::getClass<TrkrClusterContainer>(topNode,"TRKR_CLUSTER");
  if(!m_clusters)
    {
      std::cout << PHWHERE << " Can't find cluster container, can't continue."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
  if(!m_tGeometry)
    {
      std::cout << PHWHERE << " Can't find ActsGeometry, can't continue."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}
