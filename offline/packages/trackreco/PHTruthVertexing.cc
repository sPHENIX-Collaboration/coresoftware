#include "PHTruthVertexing.h"

#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertex_v1.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// gsl
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <iostream>                            // for operator<<, basic_ostream
#include <set>                                 // for _Rb_tree_iterator, set
#include <utility>                             // for pair

class PHCompositeNode;

using namespace std;

PHTruthVertexing::PHTruthVertexing(const std::string& name)
  : PHInitVertexing(name)
  , _g4truth_container(nullptr)
  , _vertex_error({0.0005, 0.0005, 0.0005})
  , _embed_only(false)
{

  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);

  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  //  cout << Name() << " random seed: " << seed << endl;
  gsl_rng_set(m_RandomGenerator, seed);
}



PHTruthVertexing::~PHTruthVertexing() {
  gsl_rng_free(m_RandomGenerator);
}

int PHTruthVertexing::Setup(PHCompositeNode* topNode)
{
  int ret = PHInitVertexing::Setup(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthVertexing::Process(PHCompositeNode* topNode)
{
  /*
  if(Verbosity() > 1)
    {
      std::cout << "embed only: " << std::boolalpha << _embed_only << std::endl;
    }
  */

  // we just copy all vertices from the truth container to the SvtxVertexMap

  auto vrange =  _g4truth_container->GetPrimaryVtxRange();
  for (auto iter = vrange.first; iter != vrange.second; ++iter)  // process all primary vertexes
    {
      const int point_id = iter->first;
      int gembed =  _g4truth_container->isEmbededVtx(point_id);

      std::vector<float> pos;
      pos.clear();
      pos.assign(3, 0.0);
            
      pos[0] = iter->second->get_x();
      pos[1] = iter->second->get_y();
      pos[2] = iter->second->get_z();

      if(Verbosity() > 1)
	{
	  cout << " PHTruthVertexing: gembed: " << gembed << " nominal position: " << " vx " << pos[0] << " vy " << pos[1] << " vz " << pos[2] << endl;
	}

      // We take all primary vertices and copy them to the vertex map, regardless of track embedding ID 

      pos[0] += _vertex_error[0] * gsl_ran_ugaussian(m_RandomGenerator);
      pos[1] += _vertex_error[1] * gsl_ran_ugaussian(m_RandomGenerator);
      pos[2] += _vertex_error[2] * gsl_ran_ugaussian(m_RandomGenerator);
      
      
      if (Verbosity() > 0)
	{
	  cout << __LINE__ << " PHTruthVertexing::Process: point_id " << point_id << " gembed " << gembed << "   {" << pos[0]
	       << ", " << pos[1] << ", " << pos[2] << "} +- {"
	       << _vertex_error[0] << ", " << _vertex_error[1] << ", "
	       << _vertex_error[2] << "}" << endl;
	}
      
      SvtxVertex* vertex = new SvtxVertex_v1();
      
      vertex->set_x(pos[0]);
      vertex->set_y(pos[1]);
      vertex->set_z(pos[2]);
      
      for (int j = 0; j < 3; ++j)
	{
	  for (int i = j; i < 3; ++i)
	    {
	      vertex->set_error(i, j,
				(i == j ? _vertex_error[i] * _vertex_error[i] : 0));
	    }
	}
      
      vertex->set_id(0);
      vertex->set_t0(0);
      vertex->set_chisq(0);
      vertex->set_ndof(1);
      _vertex_map->insert(vertex);
      
    }
  
  if (Verbosity() > 0) 
      _vertex_map->identify();
   
  if(_associate_tracks)
    { assignTracksVertices(topNode); }
    
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHTruthVertexing::assignTracksVertices(PHCompositeNode *topNode)
{
  SvtxTrackMap * trackMap = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name.c_str());
  if(!trackMap)
    {
      std::cout << PHWHERE 
		<< "Can't find requested track map. Exiting" 
		<< std::endl;
      return;
    }

  for(SvtxTrackMap::Iter iter = trackMap->begin();
      iter != trackMap->end();
      ++iter)
    {
      auto trackKey = iter->first;
      auto track = iter->second;

      if(Verbosity() > 3)
	track->identify();
    
      const double trackZ = track->get_z();
      
      double dz = 9999.;
      int trackVertexId = 9999;
      
      for(SvtxVertexMap::Iter viter = _vertex_map->begin();
	  viter != _vertex_map->end();
	  ++viter)
	{
	  auto vertexKey = viter->first;
	  auto vertex = viter->second;
	  if(Verbosity() > 3)
	    vertex->identify();
	  
	  const double vertexZ = vertex->get_z();
	  
	  if( fabs(trackZ - vertexZ) < dz )
	    {
	      dz = fabs(trackZ - vertexZ);
	      trackVertexId = vertexKey;
	    }

	}

      track->set_vertex_id(trackVertexId);
      auto vertex = _vertex_map->find(trackVertexId)->second;
      vertex->insert_track(trackKey);
      
    }


}

int PHTruthVertexing::GetNodes(PHCompositeNode* topNode)
{
  _g4truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_g4truth_container)
  {
    cerr << PHWHERE << " ERROR: Can't find node G4TruthInfo" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthVertexing::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
