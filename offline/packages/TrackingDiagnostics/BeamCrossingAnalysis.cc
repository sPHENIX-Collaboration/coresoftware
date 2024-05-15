#include "BeamCrossingAnalysis.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <trackbase/TrackVertexCrossingAssoc.h>

#include <utility>

#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>

BeamCrossingAnalysis::BeamCrossingAnalysis(const std::string& name)
  : SubsysReco(name)
{
}

int BeamCrossingAnalysis::InitRun(PHCompositeNode* topNode)
{
  const char* cfilepath = filepath.c_str();
  fout = new TFile(cfilepath, "recreate");
  ntp_vertex = new TNtuple("ntp_vertex", "Vertex ntuple","crossing:event:vtxid:x:y:z:ntracks");
  ntp_track = new TNtuple("ntp_track", "Track ntuple", "crossing:event:trackid:quality:x:y:z:px:py:pz");

  getNodes(topNode);

  hcross = new TH1D("hcross", "Crossings", 1000, -150, 250);  // root histogram arguments: name,title,bins,minvalx,maxvalx
  hvertz = new TH1D("hvertz", "Vertex z", 1000, -20, 20);  
  htrackz = new TH1D("htrackz", "Track z", 1000, -20, 20);  
  hvertcross = new TH1D("hvertcross", "Vertex crossing", 1000, -150, 250);  
  htrackcross = new TH1D("htrackcross", "Track crossing", 1000, -150, 250);  

  _event = 0;

  return 0;
}

int BeamCrossingAnalysis::process_event(PHCompositeNode* /**topNode*/)
{

  std::set<short int> crossing_set = m_track_vertex_crossing_map->getCrossings();

  for (const auto& cross : crossing_set)
    {
      if(Verbosity() > 0)
	{

	  std::cout  << "Event " << _event << " Crossing " << cross << std::endl;
	}

      // get all vertices for this crossing 
      auto vtxit = m_track_vertex_crossing_map->getVertices(cross);

      for (auto itr = vtxit.first; itr != vtxit.second; ++itr)
	{
	  unsigned int vtxid = itr->second;	  
	  const SvtxVertex* svtxVertex = m_vertexMap->get(vtxid);
	  if(!svtxVertex) { continue; }

	  float vx=0;
	  float vy=0;
	  float vz=0;
          unsigned int ntracks = 0;
	  vx = svtxVertex->get_x() ;
	  vy = svtxVertex->get_y() ;
	  vz = svtxVertex->get_z();
	  
	  ntracks = svtxVertex->size_tracks();
	  	  
	  // fill vertex ntuple
          float vertex_info[] = {(float) cross, (float) _event,  (float) vtxid, vx, vy, vz, (float) ntracks };
	  ntp_vertex->Fill(vertex_info);

	  hvertz->Fill(vz);
	  hvertcross->Fill(cross);

	  if(Verbosity() > 0)
	    {
	      std::cout << "    crossing  " << cross << " vtxid " << vtxid  << " ntracks " << ntracks
			<< " x " << vx << " y " << vy << " z " << vz << " cm " <<  std::endl;
	    }
	} 
                      
       // get all tracks for this crossing
      unsigned int tcount = 0;
      auto trit = m_track_vertex_crossing_map->getTracks(cross);
      for (auto itr = trit.first; itr != trit.second; ++itr)
	{	  
	  unsigned int  trid = itr->second;	  

	  const SvtxTrack *track = m_svtxTrackMap->get(trid);	  
          if (!track) { continue; }

	  float tx=0;
	  float ty=0;
	  float tz=0;
	  float px=0;
	  float py=0;
	  float pz=0;


	  tx = track->get_x() ;
	  ty = track->get_y() ;
	  tz = track->get_z() ;

	  px = track->get_px();
	  py = track->get_py();
	  pz = track->get_pz();

	  float quality = track->get_quality();

	  // fill track ntuple
          float track_info[] = {(float) cross, (float) _event, (float) trid, quality, tx, ty, tz, px, py, pz };
	  ntp_track->Fill(track_info);

	  if(Verbosity() > 0)
	    {
	      std::cout << "    cross  " << cross << " evt " << _event << " trid " << trid 
			<< " qual " << quality << " x " << tx << " y " << ty << " z " << tz 
			<< " px " << px << " py " << py << " pz " << pz   << std::endl;
	    }

	  // count tracks per crossing
	  tcount++;

	  htrackz->Fill(tz);
	  htrackcross->Fill(cross);
	} 

      if(Verbosity() > 0)
	{
	  std::cout << "Crossing " << cross << " tcount " << tcount << std::endl;
	}
    
      if(tcount > 0)
	{
	  hcross->Fill(cross);  
	}
  
      // code here to identify single track crossings and crossings with no vertex (if different)


    }
  
  _event++;
  
  return 0;
}


int BeamCrossingAnalysis::getNodes(PHCompositeNode* topNode)
{
  m_svtxTrackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_svtxTrackMap)
    {
      std::cout << PHWHERE << "No SvtxTrackMap on node tree, exiting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!m_vertexMap)
    {
      std::cout << PHWHERE << "No vertexMap on node tree, exiting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_track_vertex_crossing_map = findNode::getClass<TrackVertexCrossingAssoc>(topNode,"TrackVertexCrossingAssocMap");
  if(!m_track_vertex_crossing_map)
    {
      std::cout << PHWHERE << "No TrackVertexCrossingAssocMap on node tree, exiting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

int BeamCrossingAnalysis::End(PHCompositeNode* /**topNode*/)
{
  fout->cd();
  ntp_vertex->Write();
  ntp_track->Write();
  hcross->Write();
  htrackz->Write();
  hvertz->Write();
  htrackcross->Write();
  hvertcross->Write();
  fout->Close();

  return 0;
}
