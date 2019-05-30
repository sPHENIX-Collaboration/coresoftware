#include "PHG4SvtxDeadArea.h"

#include <trackbase_historic/SvtxHitMap.h>
#include <trackbase_historic/SvtxHit.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                     // for SubsysReco

#include <phool/PHRandomSeed.h>
#include <phool/PHTimeServer.h>                     // for PHTimeServer, PHT...
#include <phool/PHTimer.h>                          // for PHTimer
#include <phool/getClass.h>
#include <phool/phool.h>                            // for PHWHERE

#include <gsl/gsl_rng.h>

#include <cstddef>                                 // for NULL
#include <iostream>
#include <vector>                                   // for vector

class PHCompositeNode;

using namespace std;

PHG4SvtxDeadArea::PHG4SvtxDeadArea(const string &name) :
  SubsysReco(name),
  _hits(NULL),
  _timer(PHTimeServer::get()->insert_new(name)) 
{
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  return;
}

PHG4SvtxDeadArea::~PHG4SvtxDeadArea()
{
  gsl_rng_free (RandomGenerator);
}

int PHG4SvtxDeadArea::InitRun(PHCompositeNode* topNode) {

  // get node containing the digitized hits
  _hits = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
  if (!_hits) {
    cout << PHWHERE << "ERROR: Can't find node SvtxHitMap" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  // setup our random number generator
  unsigned int seed = PHRandomSeed(); // fixed seed handled in PHRandomSeed()
  gsl_rng_set(RandomGenerator,seed);

  GenericFillDeadAreaMap(topNode,"MVTX");
  GenericFillDeadAreaMap(topNode,"SVTX");
  GenericFillDeadAreaMap(topNode,"INTT");

  if (Verbosity() > 0) {
    cout << "====================== PHG4SvtxDeadArea::InitRun() ========================" << endl;
    cout << " Random number seed = " << seed << endl;    
    for (std::map<int,float>::iterator iter = _eff_by_layer.begin();
	 iter != _eff_by_layer.end();
	 ++iter) {
      cout << " Active area in Layer #" << iter->first << " set to " << iter->second*100.0 << "%" << endl;
    }
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxDeadArea::process_event(PHCompositeNode *topNode) {

  _timer.get()->restart();
 
  std::vector<unsigned int> remove_hits;
  
  for (SvtxHitMap::Iter iter = _hits->begin();
       iter != _hits->end();
       ++iter) {
    SvtxHit* hit = iter->second;

    if (gsl_rng_uniform_pos(RandomGenerator) > get_hit_efficiency(hit->get_layer())) {
      remove_hits.push_back(hit->get_id());
      if(Verbosity() > 5)
	cout << "removing hit" << hit->get_id() << endl;
    }
  }

  for (unsigned int i = 0; i < remove_hits.size(); ++i) {
    _hits->erase(remove_hits[i]);
  }
  
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

void
PHG4SvtxDeadArea::GenericFillDeadAreaMap(PHCompositeNode *topNode, const std::string &detectorname)
{
  string nodename = "CYLINDERGEOM_" + detectorname;
  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode,nodename);
    
  if (!geom_container) return;

  if(Verbosity() > 0)
    cout << "Found " << nodename << endl;

  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) {
    int layer = layeriter->second->get_layer();

    if (_eff_by_layer.find(layer) == _eff_by_layer.end()) {
      _eff_by_layer.insert(std::make_pair(layer,1.0));
    }    
  }
  
  return;
}

