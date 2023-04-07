#include "PHG4MvtxTruthClusterizer.h"
#include "PHG4MvtxDigitizer.h"

#include <trackbase/MvtxDefs.h>                   

#include <trackbase/TrkrHit.h>  // for TrkrHit
#include <trackbase/TrkrHitv2.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4tracking/TrkrTruthTrack.h>
#include <mvtx/CylinderGeom_Mvtx.h>
#include <g4tracking/TrkrTruthTrackContainer.h>
#include <iostream>
#include <trackbase/TrkrHitSet.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <phool/PHIODataNode.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <set>

using std::set;
using std::map;
using std::multimap;

int PHG4MvtxTruthClusterizer::clusterize_hits(TrkrClusterContainer* clusters)
{
  // turn m_hits into clusters, filled into clusters
  std::cout << " clusters " << (clusters == nullptr) << std::endl;
  _D__DigitizeMvtxLadderCells(); // adapted from PHG4MvtxDigitizer
  _P__MvtxHitPruner();           // adapted from mvtx/MvtxHitPruner
  _C__ClusterMvtx(clusters);     // adapted from mvtx/MvtxClusterizer
  return Fun4AllReturnCodes::EVENT_OK;
}


PHG4MvtxTruthClusterizer::PHG4MvtxTruthClusterizer ( )
  : _energy_threshold    ( 0.95e-6 )
{ };

void PHG4MvtxTruthClusterizer::init_run(PHCompositeNode*& _topNode, int _verbosity) {
  init_clusterizer_base(_topNode, _verbosity);

  // from PHG4MvtxDigitizer (_D_)
  _D__InitRun(_topNode);
  // nothing to do intialize for MvtxHitPruner
}

void PHG4MvtxTruthClusterizer::check_g4hit(PHG4Hit* hit) {
  if (Verbosity()>10) std::cout << " -> Checking PHG4Hit" << std::endl;
  check_g4hit_status(hit);
  if (m_was_emb) {
    if (Verbosity()>3) {
      std::cout << PHWHERE << std::endl 
      << " -> Pre clustering " << (int) m_hits->size() << " hits" << std::endl;
    }
    TrkrClusterContainerv4 clusters{};
    clusterize_hits   (&clusters);
    transfer_clusters (&clusters);
    if (Verbosity()>3) {
      std::cout << PHWHERE << std::endl 
      << " -> Clustered " << (int) clusters.size() << " clusters" << std::endl;
    }
  }
  if (m_is_new_track) update_track();
}

void PHG4MvtxTruthClusterizer::end_of_event() {
  check_g4hit(nullptr); // flush out last data if ended in truth track
  m_hitsetkey_cnt.clear();
  if (Verbosity()>2) { 
      std::cout << PHWHERE << " :: tracks with clusters after clustering in MVTX" << std::endl;
      for (auto& track : m_truthtracks->getMap()) {
        std::cout << "  track("<< track.first <<") nclusters: " << track.second->getClusters().size();
        for (auto& cluster : track.second->getClusters()) std::cout << " " << (int) TrkrDefs::getLayer(cluster);
        std::cout << std::endl;
      }
  }
}

int PHG4MvtxTruthClusterizer::_D__InitRun(PHCompositeNode *topNode)
{
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

  _D__CalculateMvtxLadderCellADCScale(topNode);

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    std::cout << "====================== PHG4MvtxTruthClusterizer copy of PHG4MvtxDigitizer::InitRun() =====================" << std::endl;
    for (auto &miter : _max_adc)
    {
      std::cout << " Max ADC in Layer #" << miter.first << " = " << miter.second << std::endl;
    }
    for (auto &miter : _energy_scale)
    {
      std::cout << " Energy per ADC in Layer #" << miter.first << " = " << 1.0e6 * miter.second << " keV" << std::endl;
    }
    std::cout << "===========================================================================" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4MvtxTruthClusterizer::_D__CalculateMvtxLadderCellADCScale(PHCompositeNode *topNode) 
{
  // defaults to 8-bit ADC, short-axis MIP placed at 1/4 dynamic range
  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");

  if (!geom_container) return;

  if (Verbosity()) std::cout << "Found CYLINDERGEOM_MVTX node" << std::endl;

  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for (PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    int layer = layeriter->second->get_layer();
    float thickness = (layeriter->second)->get_pixel_thickness();
    float pitch = (layeriter->second)->get_pixel_x();
    float length = (layeriter->second)->get_pixel_z();

    float minpath = pitch;
    if (length < minpath) minpath = length;
    if (thickness < minpath) minpath = thickness;
    float mip_e = 0.003876 * minpath;

    if (Verbosity())
      std::cout << "mip_e = " << mip_e << std::endl;

    if (_max_adc.find(layer) == _max_adc.end())
    {
      _max_adc[layer] = 255;
      _energy_scale[layer] = mip_e / 64;
    }
  }
  return;
}

void PHG4MvtxTruthClusterizer::_D__DigitizeMvtxLadderCells() {
  //----------
  // Get Nodes
  //----------
  auto trkrhitsetcontainer = m_hits;

  //=============
  // Get the TrkrHitSetContainer node
  // Digitization

  // We want all hitsets for the Mvtx
  TrkrHitSetContainer::ConstRange hitset_range = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::mvtxId);
  for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range.first;
       hitset_iter != hitset_range.second;
       ++hitset_iter)
  {
    // we have an iterator to one TrkrHitSet for the mvtx from the trkrHitSetContainer
    // get the hitset key so we can find the layer
    TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;
    int layer = TrkrDefs::getLayer(hitsetkey);
    if (Verbosity() > 1) std::cout << "PHG4MvtxDigitizer: found hitset with key: " << hitsetkey << " in layer " << layer << std::endl;

    // get all of the hits from this hitset
    TrkrHitSet *hitset = hitset_iter->second;
    TrkrHitSet::ConstRange hit_range = hitset->getHits();
    std::set<TrkrDefs::hitkey> hits_rm;
    for (TrkrHitSet::ConstIterator hit_iter = hit_range.first;
         hit_iter != hit_range.second;
         ++hit_iter)
    {
      TrkrHit *hit = hit_iter->second;

      // Convert the signal value to an ADC value and write that to the hit
      // Unsigned int adc = hit->getEnergy() / (TrkrDefs::MvtxEnergyScaleup *_energy_scale[layer]);
      if (Verbosity() > 0)
        std::cout << "    PHG4MvtxDigitizer: found hit with key: " << hit_iter->first << " and signal " << hit->getEnergy() / TrkrDefs::MvtxEnergyScaleup << " in layer " << layer << std::endl;
      // Remove the hits with energy under threshold
      bool rm_hit = false;
      if ((hit->getEnergy() / TrkrDefs::MvtxEnergyScaleup) < _energy_threshold)
      {
        if (Verbosity() > 0) std::cout << "         remove hit, below energy threshold of " << _energy_threshold << std::endl;
        rm_hit = true;
      }
      unsigned short adc = (unsigned short) (hit->getEnergy() / (TrkrDefs::MvtxEnergyScaleup * _energy_scale[layer]));
      if (adc > _max_adc[layer]) adc = _max_adc[layer];
      hit->setAdc(adc);

      if (rm_hit) hits_rm.insert(hit_iter->first);
    }

    for (const auto &key : hits_rm)
    {
      if (Verbosity() > 0) std::cout << "    PHG4MvtxDigitizer: remove hit with key: " << key << std::endl;
      hitset->removeHit(key);
    }
  }
  return;
}

int PHG4MvtxTruthClusterizer::_P__MvtxHitPruner() {
  // We want to combine all strobe values for a given hitset
  // Start by looping over all MVTX hitsets and making a map of physical sensor to hitsetkey-with-strobe
  //=============================================================================
  std::multimap<TrkrDefs::hitsetkey, TrkrDefs::hitsetkey> hitset_multimap;  // will map (bare hitset, hitset with strobe)
  std::set<TrkrDefs::hitsetkey> bare_hitset_set;  // list of all physical sensor hitsetkeys (i.e. with strobe set to zero)

  TrkrHitSetContainer::ConstRange hitsetrange =
    m_hits->getHitSets(TrkrDefs::TrkrId::mvtxId);
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
    {
      auto hitsetkey = hitsetitr->first;

      // get the hitsetkey value for strobe 0
      unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
      unsigned int stave =  MvtxDefs::getStaveId(hitsetitr->first);
      unsigned int chip =  MvtxDefs::getChipId(hitsetitr->first);
      auto bare_hitsetkey =  MvtxDefs::genHitSetKey(layer, stave, chip, 0);

      hitset_multimap.insert(std::make_pair(bare_hitsetkey, hitsetkey));
      bare_hitset_set.insert(bare_hitsetkey);

      if(Verbosity() > 0) std::cout << " found hitsetkey " << hitsetkey << " for bare_hitsetkey " << bare_hitsetkey << std::endl;
    }

  // Now consolidate all hits into the hitset with strobe 0, and delete the other hitsets
  //==============================================================
  for(auto bare_it = bare_hitset_set.begin(); bare_it != bare_hitset_set.end(); ++bare_it)
    {
      auto bare_hitsetkey = *bare_it;
      TrkrHitSet* bare_hitset = (m_hits->findOrAddHitSet(bare_hitsetkey))->second;
      if(Verbosity() > 0) std::cout << "         bare_hitset " << bare_hitsetkey << " initially has " << bare_hitset->size() << " hits " << std::endl; 

      auto bare_hitsetrange= hitset_multimap.equal_range(bare_hitsetkey);
      for(auto it = bare_hitsetrange.first; it != bare_hitsetrange.second; ++ it)
	{ 
	  auto hitsetkey = it->second;

	  int strobe = MvtxDefs::getStrobeId(hitsetkey);
	  if(strobe != 0)
	    {
	      if(Verbosity() > 0)  std::cout << "            process hitsetkey " << hitsetkey << " for bare_hitsetkey " << bare_hitsetkey << std::endl;

	      // copy all hits to the hitset with strobe 0
	      TrkrHitSet* hitset = m_hits->findHitSet(hitsetkey);		

	       if(Verbosity() > 0) 
		 std::cout << "                hitsetkey " << hitsetkey << " has strobe " << strobe << " and has " << hitset->size() << " hits,  so copy it" << std::endl;

	      TrkrHitSet::ConstRange hitrangei = hitset->getHits();
	      for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
		   hitr != hitrangei.second;
		   ++hitr)
		{
		  auto hitkey = hitr->first;
		  if(Verbosity() > 0) std::cout << "                 found hitkey " << hitkey << std::endl;		  
		  // if it is already there, leave it alone, this is a duplicate hit
		  auto tmp_hit = bare_hitset->getHit(hitkey);
		  if(tmp_hit) 
		    {
		      if(Verbosity() > 0) std::cout << "                          hitkey " << hitkey << " is already in bare hitsest, do not copy" << std::endl;
		      continue;
		    }

		  // otherwise copy the hit over 
		   if(Verbosity() > 0)  std::cout << "                          copying over hitkey " << hitkey << std::endl;
		  auto old_hit = hitr->second;
		  TrkrHit *new_hit = new TrkrHitv2();
		  new_hit->setAdc(old_hit->getAdc());
		  bare_hitset->addHitSpecificKey(hitkey, new_hit);
		}

	      // all hits are copied over to the strobe zero hitset, remove this hitset
	      m_hits->removeHitSet(hitsetkey);
	    }
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

// ---------------------------------------
// mvtx/MvtxClusterizer _C__
// ---------------------------------------
void PHG4MvtxTruthClusterizer::_C__ClusterMvtx(TrkrClusterContainer* m_clusterlist) {
  // already inherit m_hits from class
  if (Verbosity() > 0)
    std::cout << "Entering PHG4MvtxTruthClusterizer::_C__ MvtxClusterizer::ClusterMvtx " << std::endl;

  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(m_topNode, "CYLINDERGEOM_MVTX");
  if (!geom_container) return;

  //-----------
  // Clustering
  //-----------

  // loop over each MvtxHitSet object (chip)
  TrkrHitSetContainer::ConstRange hitsetrange =
    m_hits->getHitSets(TrkrDefs::TrkrId::mvtxId);
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
    { // hitsetitr    : pair(TrkrDefs::hitsetkey, TrkrHitSet>;   TrkrHitSet : map <HitKey, TrkrHit>
      TrkrHitSet *hitset = hitsetitr->second; // hitset : map <TrkrDefs::hitkey, TrkrHit>
      
      if(Verbosity() > 0)
	{ 
	  unsigned int layer  = TrkrDefs::getLayer    (hitsetitr ->first);
	  unsigned int stave  = MvtxDefs::getStaveId  (hitsetitr ->first);
	  unsigned int chip   = MvtxDefs::getChipId   (hitsetitr ->first);
	  unsigned int strobe = MvtxDefs::getStrobeId (hitsetitr ->first);
	  std::cout << "MvtxClusterizer found hitsetkey " << hitsetitr->first << " layer " << layer << " stave " << stave << " chip " << chip << " strobe " << strobe << std::endl;
     	}

      if (Verbosity() > 2)
	hitset->identify();
      
      TrkrHitSet::ConstRange hitrangei = hitset->getHits();

	  auto ckey = TrkrDefs::genClusKey(hitset->getHitSetKey(), 0); // there is only one cluster made per cluskey
	  
	  // determine the size of the cluster in phi and z
	  set<int> phibins;
	  set<int> zbins;
	  
	  // determine the cluster position...
	  double locxsum = 0.;
	  double loczsum = 0.;
	  
	  double locclusx = NAN;
	  double locclusz = NAN;
	  
	  // we need the geometry object for this layer to get the global positions
	  int layer = TrkrDefs::getLayer(ckey);
	  auto layergeom = dynamic_cast<CylinderGeom_Mvtx *>(geom_container->GetLayerGeom(layer));
	  if (!layergeom)
	    exit(1);
	  
	  const unsigned int nhits = std::distance( hitrangei.first, hitrangei.second );
	  for ( auto ihit = hitrangei.first; ihit != hitrangei.second; ++ihit)
	    {
	      // size
	      int col =  MvtxDefs::getCol( ihit->first);
	      int row = MvtxDefs::getRow(  ihit->first);
	      zbins.insert(col);
	      phibins.insert(row);
	      
	      // get local coordinates, in stave reference frame, for hit
	      auto local_coords = layergeom->get_local_coords_from_pixel(row,col);
	      
	      /*
		manually offset position along y (thickness of the sensor),
		to account for effective hit position in the sensor, resulting from diffusion.
		Effective position corresponds to 1um above the middle of the sensor
	      */
	      local_coords.SetY( 1e-4 );
	      
	      // update cluster position
	      locxsum += local_coords.X();
	      loczsum += local_coords.Z();
	      // add the association between this cluster key and this hitkey to the table
	    }  //mapiter
	  
	  // This is the local position
	  locclusx = locxsum / nhits;
	  locclusz = loczsum / nhits;
	  
	  const double pitch   = layergeom->get_pixel_x();
	  const double length  = layergeom->get_pixel_z();
	  const double phisize = phibins.size() * pitch;
	  const double zsize   = zbins.size()   * length;
	  
	  if(Verbosity() > 0) {
	    std::cout << " MvtxClusterizer: cluskey " << ckey << " layer " << layer << " rad " << layergeom->get_radius() << " phibins " << phibins.size() << " pitch " << pitch << " phisize " << phisize 
		 << " zbins " << zbins.size() << " length " << length << " zsize " << zsize 
		 << " local x " << locclusx << " local y " << locclusz
		 << std::endl;
    }
	  
    // ok force it use use cluster version v4 for now (Valgrind is not happy with application of v5)
	  /* if (m_cluster_version==4){ */
    if (m_cluster_version == 4) {
	    auto clus = std::make_unique<TrkrClusterv4>();
	    clus->setAdc(nhits);
	    clus->setLocalX(locclusx);
	    clus->setLocalY(locclusz);
	    
	    clus->setPhiSize(phibins.size());
	    clus->setZSize(zbins.size());
	    // All silicon surfaces have a 1-1 map to hitsetkey. 
	    // So set subsurface key to 0
	    clus->setSubSurfKey(0);
	    
	    if (Verbosity() > 2)
	      clus->identify();
	    
	    m_clusterlist->addClusterSpecifyKey(ckey, clus.release());
    } else {
      std::cout << PHWHERE << std::endl;
      std::cout << "Error: only cluster version 4 allowed." << std::endl;
    }  // loop over hitsets
  }
  return;
}

