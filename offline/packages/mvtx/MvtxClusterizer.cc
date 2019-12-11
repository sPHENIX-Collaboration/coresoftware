/**
 * @file mvtx/MvtxClusterizer.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of MvtxClusterizer
 */
#include "MvtxClusterizer.h"
#include "MvtxDefs.h"
#include "CylinderGeom_Mvtx.h"

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrDefs.h>                     // for hitkey, getLayer
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                     // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                           // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                         // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>                            // for PHWHERE

#include <TMatrixFfwd.h>                            // for TMatrixF
#include <TMatrixT.h>                               // for TMatrixT, operator*
#include <TMatrixTUtils.h>                          // for TMatrixTRow
#include <TVector3.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <cmath>
#include <cstdlib>                                 // for exit
#include <iostream>
#include <map>                                      // for multimap<>::iterator
#include <set>                                      // for set, set<>::iterator
#include <string>
#include <vector>                                   // for vector

using namespace boost;
using namespace std;


bool MvtxClusterizer::are_adjacent(const std::pair<TrkrDefs::hitkey, TrkrHit*> &lhs, const std::pair<TrkrDefs::hitkey, TrkrHit*> &rhs)
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

MvtxClusterizer::MvtxClusterizer(const string &name)
  : SubsysReco(name)
  , m_hits(nullptr)
  , m_clusterlist(nullptr)
  , m_clusterhitassoc(nullptr)
  , m_makeZClustering(true)
{
}

int MvtxClusterizer::InitRun(PHCompositeNode *topNode)
{
  //-----------------
  // Add Cluster Node
  //-----------------

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter_dst(dstNode);

  // Create the SVX node if required
  PHCompositeNode *svxNode = dynamic_cast<PHCompositeNode *>(iter_dst.findFirst("PHCompositeNode", "TRKR"));
  if (!svxNode)
  {
    svxNode = new PHCompositeNode("TRKR");
    dstNode->addNode(svxNode);
  }

  // Create the Cluster node if required
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

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    cout << "====================== MvtxClusterizer::InitRun() =====================" << endl;
    cout << " Z-dimension Clustering = " << boolalpha << m_makeZClustering << noboolalpha << endl;
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxClusterizer::process_event(PHCompositeNode *topNode)
{
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

  // run clustering
  ClusterMvtx(topNode);
  PrintClusters(topNode);

  // done
  return Fun4AllReturnCodes::EVENT_OK;
}

void MvtxClusterizer::ClusterMvtx(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    cout << "Entering MvtxClusterizer::ClusterMvtx " << endl;

 PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
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

    if(Verbosity() > 1) cout << "MvtxClusterizer found hitsetkey " << hitsetitr->first << endl;

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
    if (Verbosity() > 2)
      cout << "hitvec.size(): " << hitvec.size() << endl;

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
	pair<multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator,  
		      multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator>  clusrange = clusters.equal_range(clusid);
	multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator mapiter = clusrange.first;
	
	if (Verbosity() > 2)
	  cout << "Filling cluster id " << clusid << endl;
	
	// make the cluster directly in the node tree
	TrkrDefs::cluskey ckey = MvtxDefs::genClusKey(hitset->getHitSetKey(), clusid);
	TrkrClusterv1 *clus = static_cast<TrkrClusterv1 *>((m_clusterlist->findOrAddCluster(ckey))->second);
	
	// determine the size of the cluster in phi and z
	set<int> phibins;
	set<int> zbins;
	
	// determine the cluster position...
	double xsum = 0.0;
	double ysum = 0.0;
	double zsum = 0.0;
	unsigned nhits = 0;

	double clusx = NAN;
	double clusy = NAN;
	double clusz = NAN;

	// we need the geometry object for this layer to get the global positions
	int layer = TrkrDefs::getLayer(ckey);
	CylinderGeom_Mvtx *layergeom = dynamic_cast<CylinderGeom_Mvtx *>(geom_container->GetLayerGeom(layer));
	if (!layergeom)
	  exit(1);

	int chip = MvtxDefs::getChipId(ckey);
	int stave =  MvtxDefs::getStaveId(ckey); 
	
	for (mapiter = clusrange.first; mapiter != clusrange.second; ++mapiter)
	  {
	    // size
	    int col =  MvtxDefs::getCol( (mapiter->second).first);
	    int row = MvtxDefs::getRow( (mapiter->second).first);
	    zbins.insert(col);
	    phibins.insert(row);

	    int pixnum = layergeom->get_pixel_number_from_xbin_zbin(row, col);
	    //cout << "   new mvtx clusterizer: cluster key " << ckey << " layer " << layer << " chip " << chip << " stave " << stave 
	    //	 << " row " << row << " col " << col << " pixnum " << pixnum << endl;;

	    TVector3 local_coords = layergeom->get_local_coords_from_pixel(pixnum);
	    TVector3 world_coords = layergeom->get_world_from_local_coords(stave, 0, 0, chip, local_coords);
	    //cout << "   new: world coords: X " << world_coords.X() << " Y " << world_coords.Y() << " Z " << world_coords.Z() << endl;
	    
	    // find the center of the pixel in local coords
	    xsum += world_coords.X();
	    ysum += world_coords.Y();
	    zsum += world_coords.Z();

	    // add the association between this cluster key and this hitkey to the table
	    m_clusterhitassoc->addAssoc(ckey, mapiter->second.first);
	    
	    ++nhits;
	  }  //mapiter
	
	
	// This is the global position
	clusx = xsum / nhits;
	clusy = ysum / nhits;
	clusz = zsum / nhits;
	//cout << "new mvtx clusterizer: clusx " << clusx << " clusy " << clusy << " clusz " << clusz << endl;
	clus->setAdc(nhits);
				
	clus->setPosition(0, clusx);
	clus->setPosition(1, clusy);
	clus->setPosition(2, clusz);
	clus->setGlobal();
	
	double thickness = layergeom->get_pixel_thickness();
	double pitch = layergeom->get_pixel_x();
	double length = layergeom->get_pixel_z();
	double phisize = phibins.size() * pitch;
	double zsize = zbins.size() * length;
	
	double ladder_location[3] = {0.0, 0.0, 0.0};
	// returns the center of the sensor in world coordinates - used to get the ladder phi location
	layergeom->find_sensor_center(stave, 0, 0, chip, ladder_location);
	double ladderphi = atan2(ladder_location[1], ladder_location[0]);
	ladderphi += layergeom->get_stave_phi_tilt();
	
	// tilt refers to a rotation around the radial vector from the origin, and this is zero for the MVTX ladders
	float tilt = 0.0;
	
	double invsqrt12 = 1.0 / sqrt(12.0);
	
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
	
	TMatrixF ERR(3, 3);
	ERR[0][0] = pow(0.5 * thickness * invsqrt12, 2);
	ERR[0][1] = 0.0;
	ERR[0][2] = 0.0;
	ERR[1][0] = 0.0;
	ERR[1][1] = pow(0.5 * phisize * invsqrt12, 2);
	ERR[1][2] = 0.0;
	ERR[2][0] = 0.0;
	ERR[2][1] = 0.0;
	ERR[2][2] = pow(0.5 * zsize * invsqrt12, 2);
	
	TMatrixF ROT(3, 3);
	ROT[0][0] = cos(ladderphi);
	ROT[0][1] = -1.0 * sin(ladderphi);
	ROT[0][2] = 0.0;
	ROT[1][0] = sin(ladderphi);
	ROT[1][1] = cos(ladderphi);
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
	
	
	//cout << "MvtxClusterizer (x,y,z) = " << clusx << "  " << clusy << "  " << clusz << endl;
	
	if (Verbosity() > 2)
	clus->identify();
	
	// Add the hit associations to the TrkrClusterHitAssoc node
	// we need the cluster key and all associated hit keys
	/*
	for(unsigned int i=0;i<hitvec.size();i++)
	  {
	    m_clusterhitassoc->addAssoc(ckey, hitvec[i].first);
	  }
	*/
      }  // clusitr
  }    // hitsetitr
  
  if(Verbosity() > 1)
    {
      // check that the associations were written correctly
      m_clusterhitassoc->identify();
    }

  return;
}

void MvtxClusterizer::PrintClusters(PHCompositeNode *topNode)
{
  if (Verbosity() >= 1)
  {
    TrkrClusterContainer *clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    if (!clusterlist) return;

    cout << "================= Aftyer MvtxClusterizer::process_event() ====================" << endl;

    cout << " There are " << clusterlist->size() << " clusters recorded: " << endl;

    clusterlist->identify();

    cout << "===========================================================================" << endl;
  }

  return;
}
