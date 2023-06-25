#include "TruthMvtxClusterBuilder.h"

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4tracking/TrkrTruthTrack.h>
#include <g4tracking/TrkrTruthTrackContainer.h>
#include <mvtx/CylinderGeom_Mvtx.h>
#include <mvtx/MvtxHitPruner.h>
#include <phool/phool.h>  // for PHWHERE
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>  // // make iwyu happy
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitv2.h>  // for TrkrHit

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#pragma GCC diagnostic pop

#include <array>
#include <cmath>
#include <cstdlib>                                 // for exit
#include <iostream>
#include <map>                                      // for multimap<>::iterator
#include <set>                                      // for set, set<>::iterator
#include <string>
#include <vector>                                   // for vector

// using statements
//   (offline/packages/mvtx/MvtxClusterizer.cc has 
//    just `using namespace boost;` and `using namespace std;`)
using boost::add_edge;
using boost::concepts::Graph;
using boost::concepts::Graph;
using boost::connected_components;
using boost::num_vertices;
using boost::undirectedS;
using boost::adjacency_list;
/* using boost::are_adjacent; */
using boost::vecS;

using std::cout;
using std::endl;
using std::make_pair;
using std::set;
using std::map;
using std::multimap;
using std::vector;

TruthMvtxClusterBuilder::TruthMvtxClusterBuilder(TrkrClusterContainer* _clusters, 
    TrkrTruthTrackContainer* _truth_tracks,
    int _verbosity) :
    m_hits        { new TrkrHitSetContainerv1() },
    m_clusters    { _clusters },
    m_truthtracks { _truth_tracks },
    m_verbosity   { _verbosity }
{};

TruthMvtxClusterBuilder::~TruthMvtxClusterBuilder() { 
  delete m_hits; 
};

void TruthMvtxClusterBuilder::check_g4hit(PHG4Hit* hit) {
  int new_trkid = hit->get_trkid();
  if (trkid != new_trkid) { // add new id if new; return in either case
    if (is_emb) cluster_hits(); // clear out the old hits
    trkid = new_trkid;
    is_emb = m_truthinfo->isEmbeded(trkid);
  }
}

void TruthMvtxClusterBuilder::addhitset(
    TrkrDefs::hitsetkey hitsetkey, 
    TrkrDefs::hitkey hitkey, 
    float neffelectrons) 
{
  // copy of code in PHG4TpcPadPlaneReadout::MapToPadPlane, with a switch
  // to ignore non embedded tracks
  if (!is_emb) return;

  // Add the hitset to the current embedded track
  // Code from PHG4TpcPadPlaneReadout::MapToPadPlane (around lines {}.cc::386-401)
  TrkrHitSetContainer::Iterator hitsetit = m_hits->findOrAddHitSet(hitsetkey);
  // See if this hit already exists
  TrkrHit *hit = nullptr;
  hit = hitsetit->second->getHit(hitkey);
  if (!hit)
  {
    // create a new one
    hit = new TrkrHitv2();
    hitsetit->second->addHitSpecificKey(hitkey, hit);
  }
  // Either way, add the energy to it  -- adc values will be added at digitization
  hit->addEnergy(neffelectrons);
}

void TruthMvtxClusterBuilder::reset() {
  if (is_emb) {
    cluster_hits();
    is_emb = false;
    trkid = -1;
  }
  m_hits->Reset();
  hitsetkey_cnt.clear();
}

void TruthMvtxClusterBuilder::cluster_hits() {
  // logic taken from offline/packages/mvtx/MvtxClusterizer.cc starting about line 258
  if (m_verbosity > 0)
    cout << "Entering MvtxClusterizer::ClusterMvtx " << endl;

  if (m_MvtxHitPruner) {
    m_MvtxHitPruner->process_TrkrHitSetContainer(m_hits);
  }

  auto track = m_truthtracks->getTruthTrack(trkid, m_truthinfo);

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
      
      if(m_verbosity > 0)
	{ 
	  unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
	  unsigned int stave = MvtxDefs::getStaveId(hitsetitr->first);
	  unsigned int chip = MvtxDefs::getChipId(hitsetitr->first);
	  unsigned int strobe = MvtxDefs::getStrobeId(hitsetitr->first);
	  cout << "MvtxClusterizer found hitsetkey " << hitsetitr->first << " layer " << layer << " stave " << stave << " chip " << chip << " strobe " << strobe << endl;
     	}

      if (m_verbosity > 2)
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
      if (m_verbosity > 2) cout << "hitvec.size(): " << hitvec.size() << endl;

      if(m_verbosity > 0)
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
	      if (are_adjacent(hitvec[i], hitvec[j]))
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
      // loop over the componenets and make clusters
      for (set<int>::iterator clusiter = cluster_ids.begin(); clusiter != cluster_ids.end(); ++clusiter)
	{
	  int clusid = *clusiter;
	  auto clusrange = clusters.equal_range(clusid);
	  
	  if (m_verbosity > 2) cout << "Filling cluster id " << clusid << " of " << std::distance(cluster_ids.begin(),clusiter )<< endl;
	  
	  // make the cluster directly in the node tree
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
	  auto layergeom = dynamic_cast<CylinderGeom_Mvtx *>(m_geom_container->GetLayerGeom(layer));
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
	      /* m_clusterhitassoc->addAssoc(ckey, mapiter->second.first); */
	      
	    }  //mapiter
	  
	  // This is the local position
	  locclusx = locxsum / nhits;
	  locclusz = loczsum / nhits;

	  const double pitch = layergeom->get_pixel_x();
	  const double length = layergeom->get_pixel_z();
	  const double phisize = phibins.size() * pitch;
	  const double zsize = zbins.size() * length;
	  
	  /* static const double invsqrt12 = 1./std::sqrt(12); */
	  
	  // scale factors (phi direction)
	  /*
	    they corresponds to clusters of size (2,2), (2,3), (3,2) and (3,3) in phi and z
	    other clusters, which are very few and pathological, get a scale factor of 1
	    These scale factors are applied to produce cluster pulls with width unity
	  */
	  
	  /* double phierror = pitch * invsqrt12; */
	  
	  /* static constexpr std::array<double, 7> scalefactors_phi = {{ 0.36, 0.6,0.37,0.49,0.4,0.37,0.33 }}; */
	  /* if(phibins.size() == 1 && zbins.size() == 1) phierror*=scalefactors_phi[0]; */
	  /* else if(phibins.size() == 2 && zbins.size() == 1) phierror*=scalefactors_phi[1]; */
	  /* else if(phibins.size() == 1 && zbins.size() == 2) phierror*=scalefactors_phi[2]; */
	  /* else if( phibins.size() == 2 && zbins.size() == 2 ) phierror*=scalefactors_phi[0]; */
	  /* else if( phibins.size() == 2 && zbins.size() == 3 )  phierror*=scalefactors_phi[1]; */
	  /* else if( phibins.size() == 3 && zbins.size() == 2 )  phierror*=scalefactors_phi[2]; */
	  /* else if( phibins.size() == 3 && zbins.size() == 3 )  phierror*=scalefactors_phi[3]; */
	  
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
	  
	  if(m_verbosity > 0) {
	    cout << " MvtxClusterizer: cluskey " << ckey << " layer " << layer << " rad " << layergeom->get_radius() << " phibins " << phibins.size() << " pitch " << pitch << " phisize " << phisize 
		 << " zbins " << zbins.size() << " length " << length << " zsize " << zsize 
		 << " local x " << locclusx << " local y " << locclusz
		 << endl;
    }
	  
	  /* if(m_cluster_version==3){ */
	  /*   auto clus = std::make_unique<TrkrClusterv3>(); */
	  /*   clus->setAdc(nhits); */
	  /*   clus->setLocalX(locclusx); */
	  /*   clus->setLocalY(locclusz); */
	  /*   // Take the rphi and z uncertainty of the cluster */
	  /*   clus->setActsLocalError(0,0,square(phierror)); */
	  /*   clus->setActsLocalError(0,1,0.); */
	  /*   clus->setActsLocalError(1,0,0.); */
	  /*   clus->setActsLocalError(1,1,square(zerror)); */
	    
	  /*   // All silicon surfaces have a 1-1 map to hitsetkey. */ 
	  /*   // So set subsurface key to 0 */
	  /*   clus->setSubSurfKey(0); */
	    
	  /*   if (m_verbosity > 2) */
	  /*     clus->identify(); */
	    
	  /*   m_clusterlist->addClusterSpecifyKey(ckey, clus.release()); */
	  /* }else if(m_cluster_version==4){ */

    auto clus = std::make_unique<TrkrClusterv4>();
    clus->setAdc(nhits);
    clus->setLocalX(locclusx);
    clus->setLocalY(locclusz);

    clus->setPhiSize(phibins.size());
    clus->setZSize(zbins.size());
    // All silicon surfaces have a 1-1 map to hitsetkey. 
    // So set subsurface key to 0
    clus->setSubSurfKey(0);

    if (m_verbosity > 2) {
      clus->identify();
    }

    // update the key specific to this track
    auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(ckey);
    if (hitsetkey_cnt.find(hitsetkey)==hitsetkey_cnt.end()) {
      hitsetkey_cnt[hitsetkey] = 0;
    } else {
      hitsetkey_cnt[hitsetkey] +=1;
    }
    auto new_clusid = hitsetkey_cnt[hitsetkey]; 
    auto new_ckey = TrkrDefs::genClusKey(hitsetkey, new_clusid);
    m_clusters->addClusterSpecifyKey(new_ckey, clus.release());
    track->addCluster(new_ckey);
	  /* } */
	}  // clusitr loop
    }  // loop over hitsets
  if (m_verbosity>3) std::cout << "track " << track->getTrackid() << " has " << track->getClusters().size() << " mvtx clusters" << std::endl;
      
  m_hits->Reset();
  return;
};

bool TruthMvtxClusterBuilder::are_adjacent(const std::pair<TrkrDefs::hitkey, TrkrHit*> &lhs, const std::pair<TrkrDefs::hitkey, TrkrHit*> &rhs)
{
  if (GetZClustering())
  {
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
  //}

  return false;
}

void TruthMvtxClusterBuilder::print_clusters(int nclusprint) {
  cout << PHWHERE << ": content of clusters " << endl;
  auto& tmap = m_truthtracks->getMap();
  cout << " Number of tracks: " << tmap.size() << endl;
  for (auto& _pair : tmap) {
    auto& track = _pair.second;

    printf("id(%2i) phi:eta:pt(", (int)track->getTrackid());
    cout << "phi:eta:pt(";
    printf("%5.2f:%5.2f:%5.2f", track->getPhi(), track->getPseudoRapidity(), track->getPt());
      /* Form("%5.2:%5.2:%5.2", track->getPhi(), track->getPseudoRapidity(), track->getPt()) */
      //<<track->getPhi()<<":"<<track->getPseudoRapidity()<<":"<<track->getPt() 
      cout << ") nclusters(" << track->getClusters().size() <<") ";
    if (m_verbosity <= 10) { cout << endl; }
    else {
      int nclus = 0;
      for (auto cluskey : track->getClusters()) {
        cout << " " 
          << ((int) TrkrDefs::getHitSetKeyFromClusKey(cluskey)) <<":index(" <<
          ((int)  TrkrDefs::getClusIndex(cluskey)) << ")";
        ++nclus;
        if (nclusprint > 0 && nclus >= nclusprint) {
          cout << " ... "; 
          break;
        }
      }
    }
  }
  cout << PHWHERE << " ----- end of clusters " << endl;
};
