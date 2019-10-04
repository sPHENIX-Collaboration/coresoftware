#include "PHTruthVertexing.h"

#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertex_v1.h>

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
  , _vertex_error({0.01, 0.01, 0.01})
{
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

  auto vrange =  _g4truth_container->GetPrimaryVtxRange();
  set<int> gembed_set;
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

      // skip particles that are not primary
      if( sqrt(pos[0]*pos[0]+pos[1]*pos[1])  > 0.1) continue;

      // check to see if we have this vertex already      
      if(gembed_set.find(gembed) != gembed_set.end()) continue;

      gembed_set.insert(gembed);
      
      gsl_rng* RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
      unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
      //  cout << Name() << " random seed: " << seed << endl;
      gsl_rng_set(RandomGenerator, seed);
      
      pos[0] += _vertex_error[0] * gsl_ran_ugaussian(RandomGenerator);
      pos[1] += _vertex_error[1] * gsl_ran_ugaussian(RandomGenerator);
      pos[2] += _vertex_error[2] * gsl_ran_ugaussian(RandomGenerator);
      
      gsl_rng_free(RandomGenerator);
      
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

  /*
  PHG4VtxPoint* first_point = _g4truth_container->GetPrimaryVtx(
      _g4truth_container->GetPrimaryVertexIndex());

  std::vector<float> pos;

  pos.clear();
  pos.assign(3, 0.0);

  pos[0] = first_point->get_x();
  pos[1] = first_point->get_y();
  pos[2] = first_point->get_z();

  gsl_rng* RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
                                       //  cout << Name() << " random seed: " << seed << endl;
  gsl_rng_set(RandomGenerator, seed);

  pos[0] += _vertex_error[0] * gsl_ran_ugaussian(RandomGenerator);
  pos[1] += _vertex_error[1] * gsl_ran_ugaussian(RandomGenerator);
  pos[2] += _vertex_error[2] * gsl_ran_ugaussian(RandomGenerator);

  gsl_rng_free(RandomGenerator);

  if (Verbosity() > 1)
  {
    cout << __LINE__ << " PHTruthVertexing::Process: {" << pos[0]
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
  */

  return Fun4AllReturnCodes::EVENT_OK;
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
