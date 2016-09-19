/*!
        \file mSvtxEveDisplay.cxx
        \author Sookhyun Lee
        \brief reconstructed charged tracks and their clusters
        \version $Revision: 1.1 $
        \date    $Date: 07/26/2016
*/

// STL and BOOST includes
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <boost/bind.hpp>

// PHENIX includes
#include <phool/PHCompositeNode.h>
//#include <PHPoint.h>
#include <phool/getClass.h>

// ROOT and EVE includes
#include <TEveManager.h>
#include <TEveTrackPropagator.h>
#include <TEveTrack.h>
#include <TEvePointSet.h>
#include <TEveElement.h>
#include <TEveCaloData.h>
#include <TEveCalo.h>
#include <TH2F.h>

#include <g4hough/SvtxVertexMap.h>		// for maps & tpc
#include <g4hough/SvtxVertex.h>
#include <g4hough/SvtxTrackMap.h>
#include <g4hough/SvtxTrack.h>
#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxCluster.h>

#include <PHEveDisplay.h>

#include <mSvtxEveDisplay.h>

using boost::bind;

mSvtxEveDisplay::mSvtxEveDisplay(boost::shared_ptr<PHEveDisplay> dispin) :
  mPHEveModuleBase(),
  _evedisp(dispin),
  _prop(NULL),
  _svtx_tracks(NULL)     
{

  verbosity = _evedisp->get_verbosity();
  _evemanager = _evedisp->get_eve_manager();
  _prop = _evedisp->get_cnt_prop();
  _svtx_tracks = new TEveTrackList("Svtx Tracks");
  _evemanager->AddElement(_svtx_tracks,_evedisp->get_svtx_list());

}

mSvtxEveDisplay::~mSvtxEveDisplay()
{
}

void
mSvtxEveDisplay::init(PHCompositeNode* topNode)
{  
}

void 
mSvtxEveDisplay::init_run(PHCompositeNode* topNode)
{
}

bool
mSvtxEveDisplay::event(PHCompositeNode* topNode)
{
  clear();
  
  if(verbosity) std::cout<<"mSvtxEveDisplay - event() begins."<<std::endl;
  try
    {
      create_nodes(topNode);
      draw_tracks();      
    }
  catch(std::exception& e)
    {
      static bool first(true);
      if( first )
	std::cout << "mSvtxEveDisplay::event Exception: " << e.what() << std::endl;
      first = false;
    }

  return true;

}

void
mSvtxEveDisplay::end(PHCompositeNode* topNode)
{
}

void 
mSvtxEveDisplay::create_nodes(PHCompositeNode* topNode)
{
  _vertexmap = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
  if(!_vertexmap) std::cout<< "SvtxVertexMap node not found!!"<<std::endl;
  _trackmap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  if(!_trackmap) std::cout<< "SvtxTrackMap node not found!!" <<std::endl;
  _clustermap = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
  if(!_clustermap) std::cout<< "SvtxClusterMap node not found!!"<<std::endl;
  if(verbosity) std::cout<<"mSvtxEveDisplay - nodes created."<<std::endl;

}

void
mSvtxEveDisplay::draw_tracks()
{
  // svtx tracks
  for (SvtxVertexMap::Iter iter = _vertexmap->begin();
           iter != _vertexmap->end();
           ++iter) {

        SvtxVertex *vertex = iter->second;
	float vx = vertex->get_x();
	float vy = vertex->get_y();
	float vz = vertex->get_z();

	for(SvtxVertex::TrackIter iter = vertex->begin_tracks();
		iter != vertex->end_tracks();
		++iter){
	      SvtxTrack *track = _trackmap->get(*iter);
              float px = track->get_px();
              float py = track->get_py();
              float pz = track->get_pz();
	      int charge = track->get_charge();
	      //int  pid = track->get_pid();
 	      int pid = 0;//place holder
              TEveRecTrackT<double>* trk_reco = new TEveRecTrackT<double>();
              trk_reco->fV.Set(vx,
                               vy,
                               vz);
              trk_reco->fP.Set(px,
                               py,
                               pz);
              trk_reco->fSign =charge;
              if(!pid_cut(pid)) continue;
              TEveTrack* trk = new TEveTrack(trk_reco, _prop);
              trk->SetSmooth(kTRUE);
              if(charge > 0)
              trk->SetLineColor(kYellow);
              else
              trk->SetLineColor(kCyan);
              trk->SetLineWidth(2.5);
              for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
                       iter != track->end_clusters();
                       ++iter) {
          	     unsigned int cluster_id = *iter;
          	     SvtxCluster* cluster = _clustermap->get(cluster_id);
          	     //unsigned int layer = cluster->get_layer();
		     //float cluste = cluster->get_e();
		     //if(e < 0.1) continue;
		     float posx = cluster->get_position(0);
		     float posy = cluster->get_position(1);
		     float posz = cluster->get_position(2);
		     if(verbosity>2) std::cout<<"posx= "<<posx <<", posy= "<<posy<<", posz= "<<posz<<std::endl;
                     trk->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDaughter,
                                      TEveVectorD(posx, posy, posz)));
		     trk->SetRnrPoints(kTRUE);
         	     trk->SetMarkerStyle(25); 
		     trk->SetMarkerSize(1);
		     if(charge > 0)
		     trk->SetMarkerColor(kYellow);
		     else
		     trk->SetMarkerColor(kCyan);
              }///< SvtxTrack: ClusterIter
	      trk->MakeTrack();
              if(verbosity>2) std::cout<<"mSvtxEveDisplay - track made. "<<std::endl; 	     
 	      _svtx_tracks->AddElement(trk);
	}///< SvtxVertex: TrackIter
  }///< SvtxVertexMAp

}

bool 
mSvtxEveDisplay::pid_cut(int pid)
{
if(fabs(pid)==211 || fabs(pid)==321 || fabs(pid)==11 || fabs(pid)==13 || fabs(pid)==15 || fabs(pid)==17) return true;

if (pid == 0) return true;//always return true
return false;
}

void 
mSvtxEveDisplay::clear()
{
  _prop->IncDenyDestroy(); // Keep the track propagator from being destroyed
  _svtx_tracks->DestroyElements();
}
