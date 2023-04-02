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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#pragma GCC diagnostic pop

using boost::add_edge;
using boost::concepts::Graph;
using boost::concepts::Graph;
using boost::connected_components;
using boost::num_vertices;
using boost::undirectedS;
using boost::adjacency_list;
using boost::vecS;

using std::cout;
using std::endl;
using std::make_pair;
using std::set;
using std::map;
using std::multimap;
using std::vector;
using std::make_pair;

int PHG4MvtxTruthClusterizer::clusterize_hits(TrkrClusterContainer* clusters)
{
  std::cout << " clusters " << (clusters == nullptr) << std::endl;
  // digitize MVTX
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

void PHG4MvtxTruthClusterizer::end_of_event() {
  check_g4hit(nullptr); // flush out last data if ended in truth track
  m_hitsetkey_cnt.clear();
  if (m_verbosity>2) { 
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
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _D__CalculateMvtxLadderCellADCScale(topNode);

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    cout << "====================== PHG4MvtxTruthClusterizer copy of PHG4MvtxDigitizer::InitRun() =====================" << endl;
    for (auto &miter : _max_adc)
    {
      cout << " Max ADC in Layer #" << miter.first << " = " << miter.second << endl;
    }
    for (auto &miter : _energy_scale)
    {
      cout << " Energy per ADC in Layer #" << miter.first << " = " << 1.0e6 * miter.second << " keV" << endl;
    }
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4MvtxTruthClusterizer::_D__CalculateMvtxLadderCellADCScale(PHCompositeNode *topNode) 
{
  // defaults to 8-bit ADC, short-axis MIP placed at 1/4 dynamic range
  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");

  if (!geom_container) return;

  if (Verbosity()) cout << "Found CYLINDERGEOM_MVTX node" << endl;

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
      cout << "mip_e = " << mip_e << endl;

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
    // we have an itrator to one TrkrHitSet for the mvtx from the trkrHitSetContainer
    // get the hitset key so we can find the layer
    TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;
    int layer = TrkrDefs::getLayer(hitsetkey);
    if (Verbosity() > 1) cout << "PHG4MvtxDigitizer: found hitset with key: " << hitsetkey << " in layer " << layer << endl;

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
      //unsigned int adc = hit->getEnergy() / (TrkrDefs::MvtxEnergyScaleup *_energy_scale[layer]);
      if (Verbosity() > 0)
        cout << "    PHG4MvtxDigitizer: found hit with key: " << hit_iter->first << " and signal " << hit->getEnergy() / TrkrDefs::MvtxEnergyScaleup << " in layer " << layer << std::endl;
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
      if (Verbosity() > 0) cout << "    PHG4MvtxDigitizer: remove hit with key: " << key << endl;
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

      if(Verbosity() > 0) cout << " found hitsetkey " << hitsetkey << " for bare_hitsetkey " << bare_hitsetkey << endl;
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
	      if(Verbosity() > 0)  cout << "            process hitsetkey " << hitsetkey << " for bare_hitsetkey " << bare_hitsetkey << endl;

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
bool PHG4MvtxTruthClusterizer::_C__are_adjacent(const std::pair<TrkrDefs::hitkey, TrkrHit*> &lhs, const std::pair<TrkrDefs::hitkey, TrkrHit*> &rhs)
{
  if (GetZClustering())
  {
    // column is first, row is second
    if (fabs( MvtxDefs::getCol(lhs.first) - MvtxDefs::getCol(rhs.first) ) <= 1)
    {
      if (fabs( MvtxDefs::getRow(lhs.first) - MvtxDefs::getRow(rhs.first) ) <= 1)
      {
        return true;
      }
    }
  }
  else
  {
    if (fabs( MvtxDefs::getCol(lhs.first) - MvtxDefs::getCol(rhs.first) ) == 0)
    {
      if (fabs( MvtxDefs::getRow(lhs.first) - MvtxDefs::getRow(rhs.first) ) <= 1)
      {
        return true;
      }
    }
  }

  return false;
}

void PHG4MvtxTruthClusterizer::_C__ClusterMvtx(TrkrClusterContainer* m_clusterlist) {
  // already inherit m_hits from class
  if (Verbosity() > 0)
    cout << "Entering PHG4MvtxTruthClusterizer::_C__ MvtxClusterizer::ClusterMvtx " << endl;

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
    {
      TrkrHitSet *hitset = hitsetitr->second;
      
      if(Verbosity() > 0)
	{ 
	  unsigned int layer  = TrkrDefs::getLayer    (hitsetitr ->first);
	  unsigned int stave  = MvtxDefs::getStaveId  (hitsetitr ->first);
	  unsigned int chip   = MvtxDefs::getChipId   (hitsetitr ->first);
	  unsigned int strobe = MvtxDefs::getStrobeId (hitsetitr ->first);
	  cout << "MvtxClusterizer found hitsetkey " << hitsetitr->first << " layer " << layer << " stave " << stave << " chip " << chip << " strobe " << strobe << endl;
     	}

      if (Verbosity() > 2)
	hitset->identify();
      
      // fill a vector of hits to make things easier
      std::vector <std::pair< TrkrDefs::hitkey, TrkrHit*> > hitvec;
      
      TrkrHitSet::ConstRange hitrangei = hitset->getHits();
      for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
	   hitr != hitrangei.second;
	   ++hitr)
	{
	  hitvec.push_back(make_pair(hitr->first, hitr->second));
	}
      if (Verbosity() > 2) cout << "hitvec.size(): " << hitvec.size() << endl;

      if(Verbosity() > 0)
	{
	  for (unsigned int i = 0; i < hitvec.size(); i++)
	    {
	      auto hitkey = hitvec[i].first;
	      auto row = MvtxDefs::getRow(hitkey);
	      auto col = MvtxDefs::getCol(hitkey);
	      std::cout << "      hitkey " << hitkey << " row " << row << " col " << col << std::endl; 
	    }

	}
      
       // do the clustering
      typedef adjacency_list<vecS, vecS, undirectedS> Graph;
      Graph G;
      
      // loop over hits in this chip
      for (unsigned int i = 0; i < hitvec.size(); i++)
	{
	  for (unsigned int j = 0; j < hitvec.size(); j++)
	    {
	      if (_C__are_adjacent(hitvec[i], hitvec[j]))
		add_edge(i, j, G);
	    }
	}
      
      // Find the connections between the vertices of the graph (vertices are the rawhits,
      // connections are made when they are adjacent to one another)
      vector<int> component(num_vertices(G));
      
      // this is the actual clustering, performed by boost
      connected_components(G, &component[0]);
      
      // Loop over the components(hits) compiling a list of the
      // unique connected groups (ie. clusters).
      set<int> cluster_ids;  // unique components
      //multimap<int, pixel> clusters;
      multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*> >  clusters;
      for (unsigned int i = 0; i < component.size(); i++)
	{
	  cluster_ids.insert(component[i]);
	  clusters.insert(make_pair(component[i], hitvec[i]));
	}
  int total_clusters = 0;
      for (set<int>::iterator clusiter = cluster_ids.begin(); clusiter != cluster_ids.end(); ++clusiter)
	{
	  int clusid = *clusiter;
	  auto clusrange = clusters.equal_range(clusid);
	  
	  if (Verbosity() > 2) cout << "Filling cluster id " << clusid << " of " << std::distance(cluster_ids.begin(),clusiter )<< endl;
	  
    ++total_clusters;
	  auto ckey = TrkrDefs::genClusKey(hitset->getHitSetKey(), clusid);
	  
	  // determine the size of the cluster in phi and z
	  set<int> phibins;
	  set<int> zbins;
	  
	  // determine the cluster position...
	  double locxsum = 0.;
	  double loczsum = 0.;
	  const unsigned int nhits = std::distance( clusrange.first, clusrange.second );
	  
	  double locclusx = NAN;
	  double locclusz = NAN;
	  
	  // we need the geometry object for this layer to get the global positions
	  int layer = TrkrDefs::getLayer(ckey);
	  auto layergeom = dynamic_cast<CylinderGeom_Mvtx *>(geom_container->GetLayerGeom(layer));
	  if (!layergeom)
	    exit(1);
	  
	  for ( auto mapiter = clusrange.first; mapiter != clusrange.second; ++mapiter)
	    {
	      // size
	      int col =  MvtxDefs::getCol( (mapiter->second).first);
	      int row = MvtxDefs::getRow( (mapiter->second).first);
	      zbins.insert(col);
	      phibins.insert(row);
	      
	      // get local coordinates, in stae reference frame, for hit
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
	  
	  /* static const double invsqrt12 = 1./std::sqrt(12); */
	  
	  // scale factors (phi direction)
	  /*
	    they corresponds to clusters of size (2,2), (2,3), (3,2) and (3,3) in phi and z
	    other clusters, which are very few and pathological, get a scale factor of 1
	    These scale factors are applied to produce cluster pulls with width unity
	  */
	  
	  /* double phierror = pitch * invsqrt12; */
	  
	  /* static constexpr std::array<double, 7> scalefactors_phi = {{ 0.36, 0.6,0.37,0.49,0.4,0.37,0.33 }}; */
	  /* if     ( phibins.size() == 1 && zbins.size() == 1 ) phierror*=scalefactors_phi[0]; */
	  /* else if( phibins.size() == 2 && zbins.size() == 1 ) phierror*=scalefactors_phi[1]; */
	  /* else if( phibins.size() == 1 && zbins.size() == 2 ) phierror*=scalefactors_phi[2]; */
	  /* else if( phibins.size() == 2 && zbins.size() == 2 ) phierror*=scalefactors_phi[0]; */
	  /* else if( phibins.size() == 2 && zbins.size() == 3 ) phierror*=scalefactors_phi[1]; */
	  /* else if( phibins.size() == 3 && zbins.size() == 2 ) phierror*=scalefactors_phi[2]; */
	  /* else if( phibins.size() == 3 && zbins.size() == 3 ) phierror*=scalefactors_phi[3]; */
	  
	  
	  // scale factors (z direction)
	  /*
	    they corresponds to clusters of size (2,2), (2,3), (3,2) and (3,3) in z and phi
	    other clusters, which are very few and pathological, get a scale factor of 1
	  */
	  /* static constexpr std::array<double, 4> scalefactors_z = {{ 0.47, 0.48, 0.71, 0.55 }}; */
	  /* double zerror = length*invsqrt12; */
	  /* if( zbins.size() == 2 && phibins.size() == 2 ) zerror*=scalefactors_z[0]; */
	  /* else if( zbins.size() == 2 && phibins.size() == 3 )  zerror*=scalefactors_z[1]; */
	  /* else if( zbins.size() == 3 && phibins.size() == 2 )  zerror*=scalefactors_z[2]; */
	  /* else if( zbins.size() == 3 && phibins.size() == 3 )  zerror*=scalefactors_z[3]; */
	  
	  if(Verbosity() > 0)
	    cout << " MvtxClusterizer: cluskey " << ckey << " layer " << layer << " rad " << layergeom->get_radius() << " phibins " << phibins.size() << " pitch " << pitch << " phisize " << phisize 
		 << " zbins " << zbins.size() << " length " << length << " zsize " << zsize 
		 << " local x " << locclusx << " local y " << locclusz
		 << endl;
	  
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
    }
	  /* }else if(m_cluster_version==5){ */
	    /* auto clus = std::make_unique<TrkrClusterv5>(); */
	    /* clus->setAdc(nhits); */
	    /* clus->setMaxAdc(1); */
	    /* clus->setLocalX(locclusx); */
	    /* clus->setLocalY(locclusz); */
	    /* clus->setPhiError(phierror); */
	    /* clus->setZError(zerror); */
	    /* clus->setPhiSize(phibins.size()); */
	    /* clus->setZSize(zbins.size()); */
	    /* // All silicon surfaces have a 1-1 map to hitsetkey. */ 
	    /* // So set subsurface key to 0 */
	    /* clus->setSubSurfKey(0); */
	    
	    /* if (Verbosity() > 2) */
	    /*   clus->identify(); */
	    
	    /* m_clusterlist->addClusterSpecifyKey(ckey, clus.release()); */
	  /* } */
	}  // clusitr loop
    }  // loop over hitsets
      
  return;

}

