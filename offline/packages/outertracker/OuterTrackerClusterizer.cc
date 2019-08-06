#include "OuterTrackerClusterizer.h"
#include "CylinderGeomOuterTracker.h"
#include "OuterTrackerDefs.h"

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTimeServer.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TMatrixFfwd.h>                            // for TMatrixF
#include <TMatrixT.h>                               // for TMatrixT, operator*
#include <TMatrixTUtils.h>                          // for TMatrixTRow

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>


#include <cmath>
#include <iostream>
#include <memory>
#include <set>

using namespace boost;
using namespace std;

bool OuterTrackerClusterizer::ladder_are_adjacent( const std::pair<TrkrDefs::hitkey, TrkrHit*> &lhs, const std::pair<TrkrDefs::hitkey, TrkrHit*> &rhs, const int layer)
{
  if (get_z_clustering(layer))
    {
      if (fabs( OuterTrackerDefs::getCol(lhs.first) - OuterTrackerDefs::getCol(rhs.first) ) <= 1)
	{
	  if (fabs( OuterTrackerDefs::getRow(lhs.first) - OuterTrackerDefs::getRow(rhs.first) ) <= 1)
	    {
	      return true;
	    }
	}
    }
  else
    if (fabs( OuterTrackerDefs::getCol(lhs.first) - OuterTrackerDefs::getCol(rhs.first) ) == 0)
      {
	if (fabs( OuterTrackerDefs::getRow(lhs.first) - OuterTrackerDefs::getRow(rhs.first) ) <= 1)
	  {
	    return true;
	  }
      }

  return false;
}

OuterTrackerClusterizer::OuterTrackerClusterizer(const string& name,
                                 unsigned int min_layer,
                                 unsigned int max_layer)
  : SubsysReco(name)
  , m_hits(nullptr)
  , m_clusterlist(nullptr)
  , m_clusterhitassoc(nullptr)
  , _fraction_of_mip(0.5)
  , _thresholds_by_layer()
  , _make_z_clustering()
  , _make_e_weights()
  , _timer(PHTimeServer::get()->insert_new(name))
{
}

int OuterTrackerClusterizer::InitRun(PHCompositeNode* topNode)
{
  //-----------------
  // Add Cluster Node
  //-----------------

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter_dst(dstNode);

  // Create the Cluster and association nodes if required
  TrkrClusterContainer *trkrclusters = findNode::getClass<TrkrClusterContainer>(dstNode, "TRKR_CLUSTER");
  if (!trkrclusters)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode =
      dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
      {
	DetNode = new PHCompositeNode("TRKR");
	dstNode->addNode(DetNode);
      }
    
    trkrclusters = new TrkrClusterContainer();
    PHIODataNode<PHObject> *TrkrClusterContainerNode =
      new PHIODataNode<PHObject>(trkrclusters, "TRKR_CLUSTER", "PHObject");
    DetNode->addNode(TrkrClusterContainerNode);
  }

  TrkrClusterHitAssoc *clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode,"TRKR_CLUSTERHITASSOC");
  if(!clusterhitassoc)
    {
      PHNodeIterator dstiter(dstNode);
      PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
      if (!DetNode)
	{
	  DetNode = new PHCompositeNode("TRKR");
	  dstNode->addNode(DetNode);
	}

      clusterhitassoc = new TrkrClusterHitAssoc();
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(clusterhitassoc, "TRKR_CLUSTERHITASSOC", "PHObject");
      DetNode->addNode(newNode);
    }


  //---------------------
  // Calculate Thresholds
  //---------------------

  //CalculateLadderThresholds(topNode);

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    cout << "====================== OuterTrackerClusterizer::InitRun() =====================" << endl;
    cout << " Fraction of expected thickness MIP energy = " << _fraction_of_mip << endl;
    for (std::map<int, float>::iterator iter = _thresholds_by_layer.begin();
         iter != _thresholds_by_layer.end();
         ++iter)
    {
      cout << " Cluster Threshold in Layer #" << iter->first << " = " << 1.0e6 * iter->second << " keV" << endl;
    }
    for (std::map<int, bool>::iterator iter = _make_z_clustering.begin();
         iter != _make_z_clustering.end();
         ++iter)
    {
      cout << " Z-dimension Clustering in Layer #" << iter->first << " = " << boolalpha << iter->second << noboolalpha << endl;
    }
    for (std::map<int, bool>::iterator iter = _make_e_weights.begin();
         iter != _make_e_weights.end();
         ++iter)
    {
      cout << " Energy weighting clusters in Layer #" << iter->first << " = " << boolalpha << iter->second << noboolalpha << endl;
    }
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int OuterTrackerClusterizer::process_event(PHCompositeNode* topNode)
{
  _timer.get()->restart();

  // get node containing the digitized hits
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for clusters
  m_clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterlist)
  {
    cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTER." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for cluster hit associations
  m_clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!m_clusterhitassoc)
  {
    cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTERHITASSOC" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  ClusterLadderCells(topNode);
  PrintClusters(topNode);

  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

/*
void OuterTrackerClusterizer::CalculateLadderThresholds(PHCompositeNode* topNode)
{
  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_OUTERTRACKER");
  if (!geom_container) return;

  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for (PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    // index mapping
    int layer = layeriter->second->get_layer();

    // cluster threshold
    float thickness = (layeriter->second)->get_OuterRadius() - (layeriter->second)->get_InnerRadius();
    float threshold = _fraction_of_mip * 0.003876 * thickness;

  }

  return;
}
*/

void OuterTrackerClusterizer::ClusterLadderCells(PHCompositeNode* topNode)
{
  if (Verbosity() > 0)
    cout << "Entering OuterTrackerClusterizer::ClusterLadderCells " << endl;

  //----------
  // Get Nodes
  //----------

  // get the geometry node
  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_OUTERTRACKER");
  if (!geom_container) return;

  //-----------
  // Clustering
  //-----------

  // loop over the OuterTrackerHitSet objects
  TrkrHitSetContainer::ConstRange hitsetrange =
      m_hits->getHitSets(TrkrDefs::TrkrId::outertrackerId);
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {
    // Each hitset contains only hits that are clusterizable - i.e. belong to a single sensor
    TrkrHitSet *hitset = hitsetitr->second;

    if(Verbosity() > 1) cout << "OuterTrackerClusterizer found hitsetkey " << hitsetitr->first << endl;
    if (Verbosity() > 2)
      hitset->identify();

    // we have a single hitset, get the info that identifies the sensor
    int layer = TrkrDefs::getLayer(hitsetitr->first);

    // we will need the geometry object for this layer to get the global position	
    CylinderGeomOuterTracker* geom = dynamic_cast<CylinderGeomOuterTracker*>(geom_container->GetLayerGeom(layer));

    double thickness = geom->get_thickness();
    double pitch = geom->get_pixel_phi_spacing() * geom->get_inner_radius();
    double length = geom->get_pixel_z_spacing();
    double radius = (geom->get_inner_radius() + geom->get_outer_radius() ) / 2.0;

    // fill a vector of hits to make things easier - gets every hit in the hitset
    std::vector <std::pair< TrkrDefs::hitkey, TrkrHit*> > hitvec;
    TrkrHitSet::ConstRange hitrangei = hitset->getHits();
    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
         hitr != hitrangei.second;
         ++hitr)
      {
	hitvec.push_back(make_pair(hitr->first, hitr->second));
      }
    if (Verbosity() > 2)
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
	//cout << " OuterTracker clustering: add cluster number " << clusid << endl; 
	// get all hits for this cluster ID only
	pair<multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator,  
	     multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator>  clusrange = clusters.equal_range(clusid);
	multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator mapiter = clusrange.first;
	
	// make the cluster directly in the node tree
	TrkrDefs::cluskey ckey = OuterTrackerDefs::genClusKey(hitset->getHitSetKey(), clusid);
	TrkrClusterv1 *clus = static_cast<TrkrClusterv1 *>((m_clusterlist->findOrAddCluster(ckey))->second);

	if (Verbosity() > 2)
	  cout << "Filling cluster with key " << ckey << endl;
		
	// determine the size of the cluster in phi and z, useful for track fitting the cluster
	set<int> phibins;
	set<int> zbins;

	// determine the cluster position...
	double phisum = 0.0;
	//double xsum = 0.0;
	//double ysum = 0.0;
	double zsum = 0.0;
	double clus_adc = 0.0;
	unsigned nhits = 0;
	
	for (mapiter = clusrange.first; mapiter != clusrange.second; ++mapiter)
	  {
	    // mapiter->second.first  is the hit key
	    //cout << " adding hitkey " << mapiter->second.first << endl; 
	    int col =  OuterTrackerDefs::getCol( (mapiter->second).first);
	    int row = OuterTrackerDefs::getRow( (mapiter->second).first);
	    zbins.insert(col);
	    phibins.insert(row);

	    // mapiter->second.second is the hit
	    double hit_adc = (mapiter->second).second->getAdc();
	    
	    // now get the positions from the geometry
	    
	    double phi_location, z_location;
	    geom->find_pixel_center(
				    col,
				    row,
				    phi_location,
				    z_location );
	    
	    phisum += phi_location * hit_adc;
	    zsum += z_location * hit_adc;
	  
	    clus_adc += hit_adc;
	    ++nhits;

	    // add this cluster-hit association to the association map of (clusterkey,hitkey)
	    m_clusterhitassoc->addAssoc(ckey, mapiter->second.first);

	    if (Verbosity() > 2) cout << "     nhits = " << nhits << endl;
	    if (Verbosity() > 2)
	      {
		cout << "  From  geometry object: hit phi " << phi_location << " hit z " << z_location  << endl;
	      }
	  }

	float phisize = phibins.size() * pitch;
	float zsize = zbins.size() * length;
	double tilt = 0;

	double clusphi = NAN;	
	double clusx = NAN;
	double clusy = NAN;
	double clusz = NAN;

	clusphi = phisum / clus_adc;
	clusx = radius * cos(clusphi);
	clusy = radius * sin(clusphi);
	clusz = zsum / clus_adc;
		
	// Fill the cluster fields
	clus->setAdc(clus_adc);
	clus->setPosition(0, clusx);
	clus->setPosition(1, clusy);
	clus->setPosition(2, clusz);
	clus->setGlobal();

	float invsqrt12 = 1.0 / sqrt(12.0);
	
	TMatrixF DIM(3, 3);
	DIM[0][0] = pow(0.5 * thickness, 2);
	DIM[0][1] = 0.0;
	DIM[0][2] = 0.0;
	DIM[1][0] = 0.0;
	DIM[1][1] = pow(0.5 * phisize, 2);
	DIM[1][2] = 0.0;
	DIM[2][0] = 0.0;
	DIM[2][1] = 0.0;
	DIM[2][2] = pow(0.5 * zsize, 2);
	
	const float corr_factor = 1.0;  // ladder
	
	TMatrixF ERR(3, 3);
	ERR[0][0] = pow(thickness * invsqrt12 * corr_factor, 2);
	ERR[0][1] = 0.0;
	ERR[0][2] = 0.0;
	ERR[1][0] = 0.0;
	ERR[1][1] = pow(phisize * invsqrt12 * corr_factor, 2);
	ERR[1][2] = 0.0;
	ERR[2][0] = 0.0;
	ERR[2][1] = 0.0;
	ERR[2][2] = pow(zsize * invsqrt12 * corr_factor, 2);
	
	TMatrixF ROT(3, 3);
	ROT[0][0] = cos(clusphi);
	ROT[0][1] = -1.0 * sin(clusphi);
	ROT[0][2] = 0.0;
	ROT[1][0] = sin(clusphi);
	ROT[1][1] = cos(clusphi);
	ROT[1][2] = 0.0;
	ROT[2][0] = 0.0;
	ROT[2][1] = 0.0;
	ROT[2][2] = 1.0;
	
	TMatrixF TILT(3, 3);
	TILT[0][0] = 1.0;
	TILT[0][1] = 0.0;
	TILT[0][2] = 0.0;
	TILT[1][0] = 0.0;
	TILT[1][1] = cos(tilt);
	TILT[1][2] = -1.0 * sin(tilt);
	TILT[2][0] = 0.0;
	TILT[2][1] = sin(tilt);
	TILT[2][2] = cos(tilt);
	
	TMatrixF R(3, 3);
	R = ROT * TILT;
	R = ROT;
	
	TMatrixF R_T(3, 3);
	R_T.Transpose(R);
	
	TMatrixF COVAR_DIM(3, 3);
	COVAR_DIM = R * DIM * R_T;
	
	clus->setSize(0, 0, COVAR_DIM[0][0]);
	clus->setSize(0, 1, COVAR_DIM[0][1]);
	clus->setSize(0, 2, COVAR_DIM[0][2]);
	clus->setSize(1, 0, COVAR_DIM[1][0]);
	clus->setSize(1, 1, COVAR_DIM[1][1]);
	clus->setSize(1, 2, COVAR_DIM[1][2]);
	clus->setSize(2, 0, COVAR_DIM[2][0]);
	clus->setSize(2, 1, COVAR_DIM[2][1]);
	clus->setSize(2, 2, COVAR_DIM[2][2]);
	
	TMatrixF COVAR_ERR(3, 3);
	COVAR_ERR = R * ERR * R_T;
	
	clus->setError(0, 0, COVAR_ERR[0][0]);
	clus->setError(0, 1, COVAR_ERR[0][1]);
	clus->setError(0, 2, COVAR_ERR[0][2]);
	clus->setError(1, 0, COVAR_ERR[1][0]);
	clus->setError(1, 1, COVAR_ERR[1][1]);
	clus->setError(1, 2, COVAR_ERR[1][2]);
	clus->setError(2, 0, COVAR_ERR[2][0]);
	clus->setError(2, 1, COVAR_ERR[2][1]);
	clus->setError(2, 2, COVAR_ERR[2][2]);

	// Add the hit associations to the TrkrClusterHitAssoc node
	// we need the cluster key and all associated hit keys
	/*
	for(unsigned int i=0;i<hitvec.size();i++)
	  {
	    m_clusterhitassoc->addAssoc(ckey, hitvec[i].first);
	  }
	*/
	
      } // end loop over cluster ID's
  }  // end loop over hitsets


  if(Verbosity() > 1)
    {
      // check that the associations were written correctly
      cout << "After OuterTrackerClusterizer, cluster-hit associations are:" << endl;
      m_clusterhitassoc->identify();
    }
  
  return;
}

void OuterTrackerClusterizer::PrintClusters(PHCompositeNode* topNode)
{
  if (Verbosity() >= 1)
  {
    TrkrClusterContainer *clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    if (!clusterlist) return;

    cout << "================= After OuterTrackerClusterizer::process_event() ====================" << endl;

    cout << " There are " << clusterlist->size() << " clusters recorded: " << endl;

    clusterlist->identify();

    cout << "===========================================================================" << endl;
  }

  return;
}

