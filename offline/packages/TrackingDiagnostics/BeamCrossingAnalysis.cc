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
  ntp_vertex = new TNtuple("ntp_vertex", "Vertex ntuple","crossing:vtxid:x:y:z:ntracks");
  ntp_track = new TNtuple("ntp_track", "Track ntuple", "crossing:trackid:x:y:z:px:py:pz");

  getNodes(topNode);

  hcross = new TH1D("hcross", "Crossings", 1000, -20, 200);  // root histogram arguments: name,title,bins,minvalx,maxvalx

  return 0;
}

int BeamCrossingAnalysis::process_event(PHCompositeNode* /**topNode*/)
{
  std::set<short int> crossing_set = m_track_vertex_crossing_map->getCrossings();

  for (const auto& cross : crossing_set)
    { 
      std::cout  << "Crossing " << cross << std::endl;

      hcross->Fill(cross);  

      // get all vertices for this crossing 
      auto vtxit = m_track_vertex_crossing_map->getVertices(cross);

      for (auto itr = vtxit.first; itr != vtxit.second; ++itr)
	{
	  unsigned int vtxid = itr->second;	  
	  const SvtxVertex* svtxVertex = m_vertexMap->get(vtxid);
float vx=0;
float vy=0;
float vz=0;
          unsigned int ntracks = 0;
	  if (svtxVertex)
	    {
	      vx = svtxVertex->get_x() * Acts::UnitConstants::cm;
	      vy = svtxVertex->get_y() * Acts::UnitConstants::cm;
	      vz = svtxVertex->get_z() * Acts::UnitConstants::cm;

	      ntracks = svtxVertex->size_tracks();
	    }	 
	  
	  // fill vertex ntuple
	  float vertex_info[] = {(float) cross, (float) vtxid, vx, vy, vz, (float) ntracks };
	  ntp_vertex->Fill(vertex_info);

	  std::cout << "    crossing  " << cross << " vertex ID " << vtxid << " x " << vx << " y " << vy << " z " << vz << std::endl;
	} 
      
                
      // get all tracks for this crossing
      auto trit = m_track_vertex_crossing_map->getTracks(cross);
      for (auto itr = trit.first; itr != trit.second; ++itr)
	{	  
	  unsigned int  trid = itr->second;	  

	  const SvtxTrack *track = m_svtxTrackMap->get(trid);	  
          if (!track) continue;

float tx=0;
float ty=0;
float tz=0;
float px=0;
float py=0;
float pz=0;


           tx = track->get_x() * Acts::UnitConstants::cm;
	   ty = track->get_y() * Acts::UnitConstants::cm;
	   tz = track->get_z() * Acts::UnitConstants::cm;

	   px = track->get_px();
	   py = track->get_py();
	   pz = track->get_pz();

	  // fill track ntuple
	  float track_info[] = {(float) cross, (float) trid, tx, ty, tz, px, py, pz };
	  ntp_track->Fill(track_info);

	  std::cout << "    crossing  " << cross << " track ID " << trid << " x " << tx << " y " << ty << " z " << tz 
		    << " px " << px << " py " << py << " pz " << pz   << std::endl;

	} 
      
    }
  
  
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
  fout->Close();

  return 0;
}
