#include "PHG4InttTruthClusterizer.h"
#include "PHG4InttDigitizer.h"
#include "InttDeadMap.h"

#include <intt/CylinderGeomIntt.h>
#include <fun4all/Fun4AllBase.h>  // for Fun4AllBase::VERBOSITY_MORE
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4tracking/TrkrTruthTrack.h>
#include <g4tracking/TrkrTruthTrackContainer.h>
#include <intt/InttClusterizer.h>
#include <iostream>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <trackbase/InttDefs.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrClusterv5.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>  // for TrkrHit
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <TSystem.h>
#include <cassert>
#include <cfloat>
#include <cstdlib>  // for exit
#include <iostream>
#include <memory>  // for allocator_traits<...
#include <set>
#include <type_traits>  // for __decay_and_strip...

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#pragma GCC diagnostic pop

using namespace boost;
#include <boost/graph/connected_components.hpp>

#include <array>
#include <cmath>
#include <iostream>
#include <set>
#include <vector> 


int PHG4InttTruthClusterizer::clusterize_hits(TrkrClusterContainer* clusters)
{
  _D__DigitizeLadderCells(m_topNode); // adapted from g4intt/PHG4InttDigitizer
  _C__ClusterLadderCells(m_topNode, clusters);  // adapted from intt/InttClusterizer

  return Fun4AllReturnCodes::EVENT_OK;
}

PHG4InttTruthClusterizer::PHG4InttTruthClusterizer ( ) { };

void PHG4InttTruthClusterizer::init_run(PHCompositeNode*& _topNode, int _verbosity) {
  init_clusterizer_base(_topNode, _verbosity);
  // from PHG4MvtxDigitizer (_D_)
  _D__InitRun(_topNode);
  // nothing to do intialize for MvtxHitPruner
}

void PHG4InttTruthClusterizer::check_g4hit(PHG4Hit* hit) {
  if (m_verbosity>10) std::cout << " -> Checking PHG4Hit" << std::endl;
  check_g4hit_status(hit);
  if (m_was_emb) {
    if (m_verbosity>3) {
      std::cout << PHWHERE << std::endl 
      << " -> Pre clustering " << (int) m_hits->size() << " hits" << std::endl;
    }
    TrkrClusterContainerv4 clusters{};
    clusterize_hits   (&clusters);
    transfer_clusters (&clusters);
    if (m_verbosity>3) {
      std::cout << PHWHERE << std::endl 
      << " -> Clustered " << (int) clusters.size() << " clusters" << std::endl;
    }
  }
  if (m_is_new_track) update_track();
}

void PHG4InttTruthClusterizer::end_of_event() {
  check_g4hit(nullptr); // flush out last data if ended in truth track
  m_hitsetkey_cnt.clear();
  if (m_verbosity>2) { 
      std::cout << PHWHERE << " :: tracks with clusters after clustering in Intt" << std::endl;
      for (auto& track : m_truthtracks->getMap()) {
        std::cout << "  track("<< track.first <<") nclusters: " << track.second->getClusters().size();
        for (auto& cluster : track.second->getClusters()) std::cout << " " << (int) TrkrDefs::getLayer(cluster);
        std::cout << std::endl;
      }
  }
}


  // ---------------------------------------
  // Implementation of:
  // (1) g4intt/PHG4InttDigitizer _D__
  //     note that no noise is added during the digitization
  // ---------------------------------------
int PHG4InttTruthClusterizer::_D__InitRun(PHCompositeNode *topNode)
{
  std::cout << "PHG4InttTruthClusterizer::_D__InitRun: detector = " << detector << std::endl;
  //-------------
  // Add Hit Node
  //-------------
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _D__CalculateLadderCellADCScale(topNode);

  /* Create the run and par nodes
   * These parameters are not currently used in the PHG4InttDigitizer
    if (false) { // this will be saved when PHG4InttDigitizer runs -- skip in PHG4InttTruthClusterizer
      PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
      PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));

      std::string paramnodename = "G4CELLPARAM_" + detector;
      std::string geonodename = "G4CELLGEO_" + detector;

      UpdateParametersWithMacro(); // not used 
    // save this to the run wise tree to store on DST
      PHNodeIterator runIter(runNode);
      PHCompositeNode *RunDetNode = dynamic_cast<PHCompositeNode *>(runIter.findFirst("PHCompositeNode", detector));
      if (!RunDetNode)
      {
        RunDetNode = new PHCompositeNode(detector);
        runNode->addNode(RunDetNode);
      }
      SaveToNodeTree(RunDetNode, paramnodename);
      // save this to the parNode for use
      PHNodeIterator parIter(parNode);
      PHCompositeNode *ParDetNode = dynamic_cast<PHCompositeNode *>(parIter.findFirst("PHCompositeNode", detector));
      if (!ParDetNode)
      {
        ParDetNode = new PHCompositeNode(detector);
        parNode->addNode(ParDetNode);
      }
      PutOnParNode(ParDetNode, geonodename);
    }
  
    mNoiseMean = get_double_param("NoiseMean");
    mNoiseSigma = get_double_param("NoiseSigma");
    mEnergyPerPair = get_double_param("EnergyPerPair");
  */

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    std::cout << "====================== PHG4InttTruthClusterizer::_D__InitRun() =====================" << std::endl;
    for (std::map<int, unsigned int>::iterator iter1 = _max_adc.begin();
         iter1 != _max_adc.end();
         ++iter1)
    {
      std::cout << " Max ADC in Layer #" << iter1->first << " = " << iter1->second << std::endl;
    }
    for (std::map<int, float>::iterator iter2 = _energy_scale.begin();
         iter2 != _energy_scale.end();
         ++iter2)
    {
      std::cout << " Energy per ADC in Layer #" << iter2->first << " = " << 1.0e6 * iter2->second << " keV" << std::endl;
    }
    std::cout << "===========================================================================" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4InttTruthClusterizer::_D__CalculateLadderCellADCScale(PHCompositeNode* topNode)
{
  // FPHX 3-bit ADC, thresholds are set in "set_fphx_adc_scale".

  //PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, "G4CELL_INTT");
  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");

  //if (!geom_container || !cells) return;
  if (!geom_container) return;

  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for (PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    int layer = layeriter->second->get_layer();
    if (_max_fphx_adc.find(layer) == _max_fphx_adc.end())
    {
      std::cout << "Error: _max_fphx_adc is not available." << std::endl;
      gSystem->Exit(1);
    }
    float thickness = (layeriter->second)->get_thickness();  // cm
    float mip_e = 0.003876 * thickness;                      // GeV
    _energy_scale.insert(std::make_pair(layer, mip_e));
  }
  return;
}


void PHG4InttTruthClusterizer::_D__DigitizeLadderCells(PHCompositeNode* topNode) {
  //---------------------------
  // Get common Nodes
  //---------------------------
  const InttDeadMap *deadmap = findNode::getClass<InttDeadMap>(topNode, "DEADMAP_INTT");
  if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE)
  {
    if (deadmap)
    {
      std::cout << "PHG4InttDigitizer::DigitizeLadderCells - Use deadmap ";
      deadmap->identify();
    }
    else
    {
      std::cout << "PHG4InttDigitizer::DigitizeLadderCells - Can not find deadmap, all channels enabled " << std::endl;
    }
  }

  // Get the TrkrHitSetContainer node
  auto& trkrhitsetcontainer = m_hits; // alias local to name used in PHG4InttDigitizer
  //-------------
  // Digitization
  //-------------

  // We want all hitsets for the Intt
  TrkrHitSetContainer::ConstRange hitset_range = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::inttId);
  for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range.first;
       hitset_iter != hitset_range.second;
       ++hitset_iter)
  {
    // we have an itrator to one TrkrHitSet for the intt from the trkrHitSetContainer
    // get the hitset key so we can find the layer
    TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;
    const int layer = TrkrDefs::getLayer(hitsetkey);
    const int ladder_phi = InttDefs::getLadderPhiId(hitsetkey);
    const int ladder_z = InttDefs::getLadderZId(hitsetkey);

    if (Verbosity() > 1)
    {
      std::cout << "PHG4InttDigitizer: found hitset with key: " << hitsetkey << " in layer " << layer << std::endl;
    }
    // get all of the hits from this hitset
    TrkrHitSet *hitset = hitset_iter->second;
    TrkrHitSet::ConstRange hit_range = hitset->getHits();
    std::set<TrkrDefs::hitkey> dead_hits;  // hits on dead channel
    for (TrkrHitSet::ConstIterator hit_iter = hit_range.first;
         hit_iter != hit_range.second;
         ++hit_iter)
    {
      // ++m_nCells; // not really used by PHG4InttDigitizer

      TrkrHit *hit = hit_iter->second;
      TrkrDefs::hitkey hitkey = hit_iter->first;
      int strip_col = InttDefs::getCol(hitkey);  // strip z index
      int strip_row = InttDefs::getRow(hitkey);  // strip phi index

      // Apply deadmap here if desired
      if (deadmap)
      {
        if (deadmap->isDeadChannelIntt(
                layer,
                ladder_phi,
                ladder_z,
                strip_col,
                strip_row))
        {
          // ++m_nDeadCells;     // not really used by PHG4InttDigitizer

          if (Verbosity() >= Fun4AllBase::VERBOSITY_MORE)
          {
            std::cout << "PHG4InttDigitizer::DigitizeLadderCells - dead strip at layer " << layer << ": ";
            hit->identify();
          }

          dead_hits.insert(hit_iter->first);  // store hitkey of dead channels to be remove later
          continue;
        }
      }  //    if (deadmap)

      if (_energy_scale.count(layer) > 1)
      {
        assert(!"Error: _energy_scale has two or more keys.");
      }
      const float mip_e = _energy_scale[layer];

      std::vector<std::pair<double, double> > vadcrange = _max_fphx_adc[layer];

      int adc = 0;
      for (unsigned int irange = 0; irange < vadcrange.size(); ++irange)
      {
        if (hit->getEnergy() / TrkrDefs::InttEnergyScaleup >= vadcrange[irange].first * (double) mip_e && hit->getEnergy() / TrkrDefs::InttEnergyScaleup < vadcrange[irange].second * (double) mip_e)
        {
          adc = (unsigned short) irange;
        }
      }
      hit->setAdc(adc);

      if (Verbosity() > 2)
      {
        std::cout << "PHG4InttDigitizer: found hit with layer " << layer << " ladder_z " << ladder_z << " ladder_phi " << ladder_phi
                  << " strip_col " << strip_col << " strip_row " << strip_row << " adc " << hit->getAdc() << std::endl;
      }
    }  // end loop over hits in this hitset

    // remove hits on dead channel in TRKR_HITSET and TRKR_HITTRUTHASSOC
    for (const auto &key : dead_hits)
    {
      if (Verbosity() > 2)
      {
        std::cout << " PHG4InttDigitizer: remove hit with key: " << key << std::endl;
      }
      hitset->removeHit(key);
    }
  }  // end loop over hitsets
  return;
}

void PHG4InttTruthClusterizer::set_adc_scale(const int &layer, const std::vector<double> &userrange)
{
  if (userrange.size() != nadcbins)
  {
    std::cout << "Error: vector in set_fphx_adc_scale(vector) must have eight elements." << std::endl;
    gSystem->Exit(1);
  }
  //sort(userrange.begin(), userrange.end()); // TODO, causes GLIBC error

  std::vector<std::pair<double, double> > vadcrange;
  for (unsigned int irange = 0; irange < userrange.size(); ++irange)
  {
    if (irange == userrange.size() - 1)
    {
      vadcrange.push_back(std::make_pair(userrange[irange], FLT_MAX));
    }
    else
    {
      vadcrange.push_back(std::make_pair(userrange[irange], userrange[irange + 1]));
    }
  }
  _max_fphx_adc.insert(std::make_pair(layer, vadcrange));
}



void PHG4InttTruthClusterizer::_C__ClusterLadderCells(PHCompositeNode* topNode, TrkrClusterContainer* m_clusterlist)
{
  if (Verbosity() > 0)
    std::cout << "Entering InttClusterizer::ClusterLadderCells " << std::endl;

  //----------
  // Get Nodes
  //----------

  // get the geometry node
  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if (!geom_container) return;

  //-----------
  // Clustering
  //-----------

  // loop over the InttHitSet objects
  TrkrHitSetContainer::ConstRange hitsetrange =
      m_hits->getHitSets(TrkrDefs::TrkrId::inttId); // from TruthClusterizerBase
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {
    // Each hitset contains only hits that are clusterizable - i.e. belong to a single sensor
    TrkrHitSet *hitset = hitsetitr->second;

    if(Verbosity() > 1) std::cout << "InttClusterizer found hitsetkey " << hitsetitr->first << std::endl;
    if (Verbosity() > 2)
      hitset->identify();

    // we have a single hitset, get the info that identifies the sensor
    int layer = TrkrDefs::getLayer(hitsetitr->first);
    int ladder_z_index = InttDefs::getLadderZId(hitsetitr->first);
   
    // we will need the geometry object for this layer to get the global position  
    CylinderGeomIntt* geom = dynamic_cast<CylinderGeomIntt*>(geom_container->GetLayerGeom(layer));
    float pitch = geom->get_strip_y_spacing();
    float length = geom->get_strip_z_spacing();
    
    // fill a vector of hits to make things easier - gets every hit in the hitset
    std::vector <std::pair< TrkrDefs::hitkey, TrkrHit*> > hitvec;
    TrkrHitSet::ConstRange hitrangei = hitset->getHits();
    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
         hitr != hitrangei.second;
         ++hitr)
      {
  hitvec.push_back(std::make_pair(hitr->first, hitr->second));
      }
    if (Verbosity() > 2)
      std::cout << "hitvec.size(): " << hitvec.size() << std::endl;
    
    typedef adjacency_list<vecS, vecS, undirectedS> Graph;
    Graph G;
    
    // Find adjacent strips
    for (unsigned int i = 0; i < hitvec.size(); i++)
      {
 for (unsigned int j = i + 1; j < hitvec.size(); j++)
  {
  if (_C__ladder_are_adjacent(hitvec[i], hitvec[j], layer))
  {
  add_edge(i, j, G);
  }
  }
  
  add_edge(i, i, G);
      }
    
    // Find the connections between the vertices of the graph (vertices are the rawhits,
    // connections are made when they are adjacent to one another)
    std::vector<int> component(num_vertices(G));
    
    // this is the actual clustering, performed by boost
    connected_components(G, &component[0]);
    
    // Loop over the components(hit cells) compiling a list of the
    // unique connected groups (ie. clusters).
    std::set<int> cluster_ids;  // unique components
    
    std::multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*> >  clusters;
    for (unsigned int i = 0; i < component.size(); i++)
      {
  cluster_ids.insert(component[i]); // one entry per unique cluster id
  clusters.insert(std::make_pair(component[i], hitvec[i]));  // multiple entries per unique cluster id
      }

    // loop over the cluster ID's and make the clusters from the connected hits
    for (std::set<int>::iterator clusiter = cluster_ids.begin(); clusiter != cluster_ids.end(); ++clusiter)
      {
  int clusid = *clusiter;
  //cout << " intt clustering: add cluster number " << clusid << std::endl; 
  // get all hits for this cluster ID only
  std::pair<std::multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator,  
  std::multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator>  clusrange = clusters.equal_range(clusid);
  std::multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator mapiter = clusrange.first;
  
  // make the cluster directly in the node tree
  TrkrDefs::cluskey ckey = TrkrDefs::genClusKey(hitset->getHitSetKey(), clusid);

  if (Verbosity() > 2)
  std::cout << "Filling cluster with key " << ckey << std::endl;

  // get the bunch crossing number from the hitsetkey
  /* short int crossing = InttDefs::getTimeBucketId(hitset->getHitSetKey()); */

  // determine the size of the cluster in phi and z, useful for track fitting the cluster
  std::set<int> phibins;
  std::set<int> zbins;

  // determine the cluster position...
  double xlocalsum = 0.0;
  double ylocalsum = 0.0;
  double zlocalsum = 0.0;
  unsigned int clus_adc = 0.0;
  unsigned int clus_maxadc = 0.0;
  unsigned nhits = 0;

  //std::cout << PHWHERE << " ckey " << ckey << ":" << std::endl;  
  for (mapiter = clusrange.first; mapiter != clusrange.second; ++mapiter)
  {
  // mapiter->second.first  is the hit key
  //cout << " adding hitkey " << mapiter->second.first << std::endl; 
  int col =  InttDefs::getCol( (mapiter->second).first);
  int row = InttDefs::getRow( (mapiter->second).first);
  zbins.insert(col);
  phibins.insert(row);

  // mapiter->second.second is the hit
  unsigned int hit_adc = (mapiter->second).second->getAdc();

  // now get the positions from the geometry
  double local_hit_location[3] = {0., 0., 0.};
  geom->find_strip_center_localcoords(ladder_z_index,
  row, col,
  local_hit_location);
  
  if (_make_e_weights[layer])
  {
    xlocalsum += local_hit_location[0] * (double) hit_adc;
    ylocalsum += local_hit_location[1] * (double) hit_adc;
    zlocalsum += local_hit_location[2] * (double) hit_adc;
  }
  else
  {
    xlocalsum += local_hit_location[0];
    ylocalsum += local_hit_location[1];
    zlocalsum += local_hit_location[2];
  }
  if(hit_adc > clus_maxadc)
    clus_maxadc = hit_adc;
  clus_adc += hit_adc;
  ++nhits;

  // add this cluster-hit association to the association map of (clusterkey,hitkey)
  if (Verbosity() > 2) std::cout << "     nhits = " << nhits << std::endl;
  if (Verbosity() > 2)
  {
  std::cout << "  From  geometry object: hit x " << local_hit_location[0] << " hit y " << local_hit_location[1] << " hit z " << local_hit_location[2] << std::endl;
  std::cout << "     nhits " << nhits << " clusx  = " << xlocalsum / nhits << " clusy " << ylocalsum / nhits << " clusz " << zlocalsum / nhits << " hit_adc " << hit_adc << std::endl;
  
  }
  }

  static const float invsqrt12 = 1./sqrt(12);
  
  // scale factors (phi direction)
  /*
  they corresponds to clusters of size 1 and 2 in phi
  other clusters, which are very few and pathological, get a scale factor of 1
  These scale factors are applied to produce cluster pulls with width unity
  */

  float phierror = pitch * invsqrt12;
  
  static constexpr std::array<double, 3> scalefactors_phi = {{ 0.85, 0.4, 0.33 }};
  if( phibins.size() == 1 && layer < 5) phierror*=scalefactors_phi[0];
  else if( phibins.size() == 2 && layer < 5) phierror*=scalefactors_phi[1];
  else if( phibins.size() == 2 && layer > 4) phierror*=scalefactors_phi[2]; 
  // z error. All clusters have a z-size of 1.
  const float zerror = length * invsqrt12;

  double cluslocaly = NAN;
  double cluslocalz = NAN;

  if (_make_e_weights[layer])
    {
      cluslocaly = ylocalsum / (double) clus_adc;
      cluslocalz = zlocalsum / (double) clus_adc;
    }
  else
    {
      cluslocaly = ylocalsum / nhits;
      cluslocalz = zlocalsum / nhits;
    }
  if(m_cluster_version==4){
    auto clus = std::make_unique<TrkrClusterv4>();
    // Fill the cluster fields
    clus->setAdc(clus_adc);
    clus->setPhiSize(phibins.size());
    clus->setZSize(1);

    if(Verbosity() > 10) clus->identify();
    
    clus->setLocalX(cluslocaly);
    clus->setLocalY(cluslocalz);
    // silicon has a 1-1 map between hitsetkey and surfaces. So set to 
    // 0
    clus->setSubSurfKey(0);
    m_clusterlist->addClusterSpecifyKey(ckey, clus.release());
    
  }else if(m_cluster_version==5){
    auto clus = std::make_unique<TrkrClusterv5>();
    clus->setAdc(clus_adc);
    clus->setMaxAdc(clus_maxadc);
    clus->setLocalX(cluslocaly);
    clus->setLocalY(cluslocalz);
    clus->setPhiError(phierror);
    clus->setZError(zerror);
    clus->setPhiSize(phibins.size());
    clus->setZSize(1);
    // All silicon surfaces have a 1-1 map to hitsetkey. 
    // So set subsurface key to 0
    clus->setSubSurfKey(0);
    
    if (Verbosity() > 2)
      clus->identify();

    m_clusterlist->addClusterSpecifyKey(ckey, clus.release());
  }
      } // end loop over cluster ID's
  }  // end loop over hitsets

  return;
}

bool PHG4InttTruthClusterizer::_C__ladder_are_adjacent( const std::pair<TrkrDefs::hitkey, TrkrHit*> &lhs, const std::pair<TrkrDefs::hitkey, TrkrHit*> &rhs, const int layer)
{
  if (get_z_clustering(layer))
    {
      if (fabs( InttDefs::getCol(lhs.first) - InttDefs::getCol(rhs.first) ) <= 1)
	{
	  if (fabs( InttDefs::getRow(lhs.first) - InttDefs::getRow(rhs.first) ) <= 1)
	    {
	      return true;
	    }
	}
    }
  else
    if (fabs( InttDefs::getCol(lhs.first) - InttDefs::getCol(rhs.first) ) == 0)
      {
	if (fabs( InttDefs::getRow(lhs.first) - InttDefs::getRow(rhs.first) ) <= 1)
	  {
	    return true;
	  }
      }

  return false;
}
