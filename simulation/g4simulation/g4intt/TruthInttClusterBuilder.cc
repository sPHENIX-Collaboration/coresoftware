#include "TruthInttClusterBuilder.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4tracking/TrkrTruthTrack.h>
#include <g4tracking/TrkrTruthTrackContainer.h>
#include <intt/CylinderGeomIntt.h>

#include <trackbase/InttDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitv2.h>  // for TrkrHit

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                         // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TMatrixFfwd.h>                            // for TMatrixF
#include <TMatrixT.h>                               // for TMatrixT, operator*
#include <TMatrixTUtils.h>                          // for TMatrixTRow


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#pragma GCC diagnostic pop


#include <array>
#include <cmath>
#include <iostream>
#include <set>
#include <vector> 

// using statements
//   (offline/packages/intt/InttClusterizer.cc has 
//    just `using namespace boost;` and `using namespace std;`)
using std::cout;
using std::endl;
using std::make_pair;
using std::multimap;
using std::pair;
using std::set;
using std::vector;
using std::vector;

using boost::add_edge;
using boost::concepts::Graph;
using boost::concepts::Graph;
using boost::connected_components;
using boost::num_vertices;
using boost::undirectedS;
using boost::adjacency_list;
using boost::vecS;

TruthInttClusterBuilder::TruthInttClusterBuilder(TrkrClusterContainer* _clusters,
      TrkrTruthTrackContainer* _truth_tracks,
      int _verbosity ) :
    m_hits        { new TrkrHitSetContainerv1() },
    m_clusters    { _clusters     },
    m_truthtracks { _truth_tracks },
    m_verbosity   { _verbosity    }
{
};

TruthInttClusterBuilder::~TruthInttClusterBuilder() {
  delete m_hits;
};

void TruthInttClusterBuilder::check_g4hit(PHG4Hit* hit) {
  int new_trkid = hit->get_trkid();
  if (trkid != new_trkid) { // add new id if new; return in either case
    if (is_emb) cluster_hits(); // clear out the old hits
    trkid = new_trkid;
    is_emb = m_truthinfo->isEmbeded(trkid);
  }
}

void TruthInttClusterBuilder::addhitset(
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

void TruthInttClusterBuilder::reset() {
  if (is_emb) {
    cluster_hits();
    is_emb = false;
    trkid = -1;
  }
  m_hits->Reset();
  hitsetkey_cnt.clear();
}

void TruthInttClusterBuilder::cluster_hits() {
  // Getting code largely from coresoftware/offline/packages/intt/InttClusterizer.cc lines 346 on
  TrkrHitSetContainer::ConstRange hitsetrange =
      m_hits->getHitSets(TrkrDefs::TrkrId::inttId);

  if (hitsetrange.first == hitsetrange.second) {
    // no hits to clusters
    return;
  }

  // Get the TruthTrack
  auto track = m_truthtracks->getTruthTrack(trkid, m_truthinfo);

  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {
    // Each hitset contains only hits that are clusterizable - i.e. belong to a single sensor
    TrkrHitSet *hitset = hitsetitr->second;

    if(m_verbosity > 1) cout << "InttClusterizer found hitsetkey " << hitsetitr->first << endl;
    if (m_verbosity > 2) hitset->identify();

    // we have a single hitset, get the info that identifies the sensor
    int layer          = TrkrDefs::getLayer(hitsetitr->first);
    int ladder_z_index = InttDefs::getLadderZId(hitsetitr->first);
   
    // we will need the geometry object for this layer to get the global position	
    CylinderGeomIntt* geom = dynamic_cast<CylinderGeomIntt*>(m_geom_container->GetLayerGeom(layer));
    /* float pitch  = geom->get_strip_y_spacing(); */
    /* float length = geom->get_strip_z_spacing(); */
    
    // fill a vector of hits to make things easier - gets every hit in the hitset
    std::vector <std::pair< TrkrDefs::hitkey, TrkrHit*> > hitvec;
    TrkrHitSet::ConstRange hitrangei = hitset->getHits();
    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
         hitr != hitrangei.second;
         ++hitr)
      {
	hitvec.push_back(make_pair(hitr->first, hitr->second));
      }
    if (m_verbosity > 2)
      cout << "hitvec.size(): " << hitvec.size() << endl;
    
    typedef adjacency_list<vecS, vecS, undirectedS> Graph;
    Graph G;
    
    // Find adjacent strips
    for (unsigned int i = 0; i < hitvec.size(); i++)
      {
	for (unsigned int j = i + 1; j < hitvec.size(); j++)
	  {
	    if (ladder_are_adjacent(hitvec[i], hitvec[j], layer))
	      {
		add_edge(i, j, G);
	      }
	  }
	
	add_edge(i, i, G);
      }
    
    // Find the connections between the vertices of the graph (vertices are the rawhits,
    // connections are made when they are adjacent to one another)
    vector<int> component(num_vertices(G));
    
    // this is the actual clustering, performed by boost
    connected_components(G, &component[0]);
    
    // Loop over the components(hit cells) compiling a list of the
    // unique connected groups (ie. clusters).
    set<int> cluster_ids;  // unique components
    
    multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*> >  clusters;
    for (unsigned int i = 0; i < component.size(); i++)
      {
	cluster_ids.insert(component[i]); // one entry per unique cluster id
	clusters.insert(make_pair(component[i], hitvec[i]));  // multiple entries per unique cluster id
      }

    // loop over the cluster ID's and make the clusters from the connected hits
    for (set<int>::iterator clusiter = cluster_ids.begin(); clusiter != cluster_ids.end(); ++clusiter)
      {
	int clusid = *clusiter;
	//cout << " intt clustering: add cluster number " << clusid << endl; 
	// get all hits for this cluster ID only
	pair<multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator,  
	     multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator>  clusrange = clusters.equal_range(clusid);
	multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator mapiter = clusrange.first;
	
	// make the cluster directly in the node tree
	TrkrDefs::cluskey ckey = TrkrDefs::genClusKey(hitset->getHitSetKey(), clusid);

	if (m_verbosity > 2)
	  cout << "Filling cluster with key " << ckey << endl;

	// get the bunch crossing number from the hitsetkey
	/* short int crossing = InttDefs::getTimeBucketId(hitset->getHitSetKey()); */

	// determine the size of the cluster in phi and z, useful for track fitting the cluster
	set<int> phibins;
	set<int> zbins;

	// determine the cluster position...
	double xlocalsum = 0.0;
	double ylocalsum = 0.0;
	double zlocalsum = 0.0;
	unsigned int clus_adc = 0.0;
	unsigned nhits = 0;

	//std::cout << PHWHERE << " ckey " << ckey << ":" << std::endl;	
	for (mapiter = clusrange.first; mapiter != clusrange.second; ++mapiter)
	  {
	    // mapiter->second.first  is the hit key
	    //cout << " adding hitkey " << mapiter->second.first << endl; 
	    int col =  InttDefs::getCol( (mapiter->second).first);
	    int row = InttDefs::getRow( (mapiter->second).first);
	    zbins.insert(col);
	    phibins.insert(row);

	    // mapiter->second.second is the hit
	    unsigned int hit_adc = (mapiter->second).second->getAdc();

	    // Add clusterkey/bunch crossing to mmap
	    /* m_clustercrossingassoc->addAssoc(ckey, crossing); */

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

	    clus_adc += hit_adc;
	    ++nhits;

	    // add this cluster-hit association to the association map of (clusterkey,hitkey)
	    /* m_clusterhitassoc->addAssoc(ckey, mapiter->second.first); */

	    if (m_verbosity > 2) cout << "     nhits = " << nhits << endl;
	    if (m_verbosity > 2)
	      {
		cout << "  From  geometry object: hit x " << local_hit_location[0] << " hit y " << local_hit_location[1] << " hit z " << local_hit_location[2] << endl;
		cout << "     nhits " << nhits << " clusx  = " << xlocalsum / nhits << " clusy " << ylocalsum / nhits << " clusz " << zlocalsum / nhits << " hit_adc " << hit_adc << endl;
		
	      }
	  }

	/* static const float invsqrt12 = 1./sqrt(12); */
	
	// scale factors (phi direction)
	/*
	  they corresponds to clusters of size 1 and 2 in phi
	  other clusters, which are very few and pathological, get a scale factor of 1
	  These scale factors are applied to produce cluster pulls with width unity
	*/

	/* float phierror = pitch * invsqrt12; */
	
	/* static constexpr std::array<double, 3> scalefactors_phi = {{ 0.85, 0.4, 0.33 }}; */
	/* if( phibins.size() == 1 && layer < 5) phierror*=scalefactors_phi[0]; */
	/* else if( phibins.size() == 2 && layer < 5) phierror*=scalefactors_phi[1]; */
	/* else if( phibins.size() == 2 && layer > 4) phierror*=scalefactors_phi[2]; */ 
	// z error. All clusters have a z-size of 1.
	/* const float zerror = length * invsqrt12; */

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
	/* if(m_cluster_version==3){ */
	/*   auto clus = std::make_unique<TrkrClusterv3>(); */
	/*   // Fill the cluster fields */
	/*   clus->setAdc(clus_adc); */
	
	/*   if(Verbosity() > 10) clus->identify(); */
	
	/*   clus->setLocalX(cluslocaly); */
	/*   clus->setLocalY(cluslocalz); */
	/*   /// silicon has a 1-1 map between hitsetkey and surfaces. So set to */
	/*     /// 0 */
	/*     clus->setSubSurfKey(0); */
	/*     clus->setActsLocalError(0,0, square(phierror)); */
	/*     clus->setActsLocalError(0,1, 0.); */
	/*     clus->setActsLocalError(1,0, 0.); */
	/*     clus->setActsLocalError(1,1, square(zerror)); */
	/*     m_clusterlist->addClusterSpecifyKey(ckey, clus.release()); */
	/* } */
	/* if(m_cluster_version==4){ */
	  auto clus = std::make_unique<TrkrClusterv4>();
	  // Fill the cluster fields
	  clus->setAdc(clus_adc);
	  clus->setPhiSize(phibins.size());
	  clus->setZSize(1);

	  if(m_verbosity > 10) clus->identify();
	  
	  clus->setLocalX(cluslocaly);
	  clus->setLocalY(cluslocalz);
	  // 0
	  clus->setSubSurfKey(0);

    // add the cluster to the trkrclustercontainer
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
     } // end loop over cluster ID's
  }  // end loop over hitsets
  // end of code from offline/packages/intt/InttClusterizer.cc

  // end cluster the hits
  m_hits->Reset();
};

bool TruthInttClusterBuilder::ladder_are_adjacent( const std::pair<TrkrDefs::hitkey, TrkrHit*> &lhs, const std::pair<TrkrDefs::hitkey, TrkrHit*> &rhs, const int layer)
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

