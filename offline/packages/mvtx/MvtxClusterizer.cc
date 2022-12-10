/**
 * @file mvtx/MvtxClusterizer.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of MvtxClusterizer
 */
#include "MvtxClusterizer.h"
#include "CylinderGeom_Mvtx.h"

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrDefs.h>                     // for hitkey, getLayer
#include <trackbase/MvtxDefs.h>                   
#include <trackbase/TrkrHitv2.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterHitAssocv3.h>

#include <trackbase/RawHit.h>
#include <trackbase/RawHitSet.h>
#include <trackbase/RawHitSetContainer.h>

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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/graph/adjacency_list.hpp>
#pragma GCC diagnostic pop

#include <boost/graph/connected_components.hpp>

#include <array>
#include <cmath>
#include <cstdlib>                                 // for exit
#include <iostream>
#include <map>                                      // for multimap<>::iterator
#include <set>                                      // for set, set<>::iterator
#include <string>
#include <vector>                                   // for vector

using namespace boost;
using namespace std;

namespace
{

  /// convenience square method
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }
}

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

bool MvtxClusterizer::are_adjacent(RawHit* lhs, RawHit* rhs)
{
  if (GetZClustering())
  {
    // column is first, row is second
    if (fabs( lhs->getPhiBin() - rhs->getPhiBin()) <= 1)//col
    {
      if (fabs( lhs->getTBin() - rhs->getTBin() ) <= 1)//Row
      {
        return true;
      }
    }
  }
  else
  {
    if (fabs(lhs->getPhiBin() - rhs->getPhiBin()  ) == 0)
    {
      if (fabs( lhs->getTBin() - rhs->getTBin() ) <= 1)
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
  , m_rawhits(nullptr)
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
  auto trkrclusters = findNode::getClass<TrkrClusterContainer>(dstNode, "TRKR_CLUSTER");
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

    trkrclusters = new TrkrClusterContainerv4;
    PHIODataNode<PHObject> *TrkrClusterContainerNode =
      new PHIODataNode<PHObject>(trkrclusters, "TRKR_CLUSTER", "PHObject");
    DetNode->addNode(TrkrClusterContainerNode);
  }

  auto clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode,"TRKR_CLUSTERHITASSOC");
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

      clusterhitassoc = new TrkrClusterHitAssocv3;
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
  if(!do_read_raw){
    m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    if (!m_hits)
      {
	cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << endl;
	return Fun4AllReturnCodes::ABORTRUN;
      }
  }else{
    // get node containing the digitized hits
    m_rawhits = findNode::getClass<RawHitSetContainer>(topNode, "TRKR_RAWHITSET");
    if (!m_rawhits)
      {
	std::cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << std::endl;
	return Fun4AllReturnCodes::ABORTRUN;
      }
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
  if(!do_read_raw){
    ClusterMvtx(topNode);
  }else{
    ClusterMvtxRaw(topNode);
  }
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
      
      if(Verbosity() > 0)
	{ 
	  unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
	  unsigned int stave = MvtxDefs::getStaveId(hitsetitr->first);
	  unsigned int chip = MvtxDefs::getChipId(hitsetitr->first);
	  unsigned int strobe = MvtxDefs::getStrobeId(hitsetitr->first);
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
      //    cout << "found cluster #: "<< clusters.size()<< endl;
      // loop over the componenets and make clusters
      for (set<int>::iterator clusiter = cluster_ids.begin(); clusiter != cluster_ids.end(); ++clusiter)
	{
	  int clusid = *clusiter;
	  auto clusrange = clusters.equal_range(clusid);
	  
	  if (Verbosity() > 2) cout << "Filling cluster id " << clusid << " of " << std::distance(cluster_ids.begin(),clusiter )<< endl;
	  
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
	      m_clusterhitassoc->addAssoc(ckey, mapiter->second.first);
	      
	    }  //mapiter
	  
	  // This is the local position
	  locclusx = locxsum / nhits;
	  locclusz = loczsum / nhits;
	  
	  const double pitch = layergeom->get_pixel_x();
	  //	std::cout << " pitch: " <<  pitch << std::endl;
	  const double length = layergeom->get_pixel_z();
	  //	std::cout << " length: " << length << std::endl;
	  const double phisize = phibins.size() * pitch;
	  const double zsize = zbins.size() * length;
	  
	  static const double invsqrt12 = 1./std::sqrt(12);
	  
	  // scale factors (phi direction)
	  /*
	    they corresponds to clusters of size (2,2), (2,3), (3,2) and (3,3) in phi and z
	    other clusters, which are very few and pathological, get a scale factor of 1
	    These scale factors are applied to produce cluster pulls with width unity
	  */
	  
	  double phierror = pitch * invsqrt12;
	  
	  static constexpr std::array<double, 7> scalefactors_phi = {{ 0.36, 0.6,0.37,0.49,0.4,0.37,0.33 }};
	  if(phibins.size() == 1 && zbins.size() == 1) phierror*=scalefactors_phi[0];
	  else if(phibins.size() == 2 && zbins.size() == 1) phierror*=scalefactors_phi[1];
	  else if(phibins.size() == 1 && zbins.size() == 2) phierror*=scalefactors_phi[2];
	  else if( phibins.size() == 2 && zbins.size() == 2 ) phierror*=scalefactors_phi[0];
	  else if( phibins.size() == 2 && zbins.size() == 3 )  phierror*=scalefactors_phi[1];
	  else if( phibins.size() == 3 && zbins.size() == 2 )  phierror*=scalefactors_phi[2];
	  else if( phibins.size() == 3 && zbins.size() == 3 )  phierror*=scalefactors_phi[3];
	  
	  
	  // scale factors (z direction)
	  /*
	    they corresponds to clusters of size (2,2), (2,3), (3,2) and (3,3) in z and phi
	    other clusters, which are very few and pathological, get a scale factor of 1
	  */
	  static constexpr std::array<double, 4> scalefactors_z = {{ 0.47, 0.48, 0.71, 0.55 }};
	  double zerror = length*invsqrt12;
	  if( zbins.size() == 2 && phibins.size() == 2 ) zerror*=scalefactors_z[0];
	  else if( zbins.size() == 2 && phibins.size() == 3 )  zerror*=scalefactors_z[1];
	  else if( zbins.size() == 3 && phibins.size() == 2 )  zerror*=scalefactors_z[2];
	  else if( zbins.size() == 3 && phibins.size() == 3 )  zerror*=scalefactors_z[3];
	  
	  if(Verbosity() > 0)
	    cout << " MvtxClusterizer: cluskey " << ckey << " layer " << layer << " rad " << layergeom->get_radius() << " phibins " << phibins.size() << " pitch " << pitch << " phisize " << phisize 
		 << " zbins " << zbins.size() << " length " << length << " zsize " << zsize 
		 << " local x " << locclusx << " local y " << locclusz
		 << endl;
	  
	  if(m_cluster_version==3){
	    auto clus = std::make_unique<TrkrClusterv3>();
	    clus->setAdc(nhits);
	    clus->setLocalX(locclusx);
	    clus->setLocalY(locclusz);
	    // Take the rphi and z uncertainty of the cluster
	    clus->setActsLocalError(0,0,square(phierror));
	    clus->setActsLocalError(0,1,0.);
	    clus->setActsLocalError(1,0,0.);
	    clus->setActsLocalError(1,1,square(zerror));
	    
	    // All silicon surfaces have a 1-1 map to hitsetkey. 
	    // So set subsurface key to 0
	    clus->setSubSurfKey(0);
	    
	    if (Verbosity() > 2)
	      clus->identify();
	    
	    m_clusterlist->addClusterSpecifyKey(ckey, clus.release());
	  }else if(m_cluster_version==4){
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
	}  // clusitr loop
    }  // loop over hitsets
      
  if(Verbosity() > 1)
    {
      // check that the associations were written correctly
      m_clusterhitassoc->identify();
    }
      
  return;
}

void MvtxClusterizer::ClusterMvtxRaw(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    cout << "Entering MvtxClusterizer::ClusterMvtx " << endl;

  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  if (!geom_container) return;

  //-----------
  // Clustering
  //-----------

  // loop over each MvtxHitSet object (chip)
  RawHitSetContainer::ConstRange hitsetrange =
    m_rawhits->getHitSets(TrkrDefs::TrkrId::mvtxId);
  for (RawHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
    {
      RawHitSet *hitset = hitsetitr->second;
      
      if(Verbosity() > 0)
	{ 
	  unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
	  unsigned int stave = MvtxDefs::getStaveId(hitsetitr->first);
	  unsigned int chip = MvtxDefs::getChipId(hitsetitr->first);
	  unsigned int strobe = MvtxDefs::getStrobeId(hitsetitr->first);
	  cout << "MvtxClusterizer found hitsetkey " << hitsetitr->first << " layer " << layer << " stave " << stave << " chip " << chip << " strobe " << strobe << endl;
     	}

      if (Verbosity() > 2)
	hitset->identify();
      
      // fill a vector of hits to make things easier
      std::vector <RawHit*> hitvec;
      
      RawHitSet::ConstRange hitrangei = hitset->getHits();
      for (RawHitSet::ConstIterator hitr = hitrangei.first;
	   hitr != hitrangei.second;
	   ++hitr)
	{
	  hitvec.push_back((*hitr));
	}
      if (Verbosity() > 2) cout << "hitvec.size(): " << hitvec.size() << endl;

      
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

      multimap<int, RawHit* >  clusters;
      for (unsigned int i = 0; i < component.size(); i++)
	{
	  cluster_ids.insert(component[i]);
	  clusters.insert(make_pair(component[i], hitvec[i]));
	}
      //    cout << "found cluster #: "<< clusters.size()<< endl;
      // loop over the componenets and make clusters
      for (set<int>::iterator clusiter = cluster_ids.begin(); clusiter != cluster_ids.end(); ++clusiter)
	{
	  int clusid = *clusiter;
	  auto clusrange = clusters.equal_range(clusid);
	  
	  if (Verbosity() > 2) cout << "Filling cluster id " << clusid << " of " << std::distance(cluster_ids.begin(),clusiter )<< endl;
	  
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
	  auto layergeom = dynamic_cast<CylinderGeom_Mvtx *>(geom_container->GetLayerGeom(layer));
	  if (!layergeom)
	    exit(1);
	  
	  for ( auto mapiter = clusrange.first; mapiter != clusrange.second; ++mapiter)
	    {
	      // size
	      int col =  (mapiter->second)->getPhiBin();
	      int row =  (mapiter->second)->getTBin();
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
	      //	      m_clusterhitassoc->addAssoc(ckey, mapiter->second.first);
	      
	    }  //mapiter
	  
	  // This is the local position
	  locclusx = locxsum / nhits;
	  locclusz = loczsum / nhits;
	  
	  const double pitch = layergeom->get_pixel_x();
	  //	std::cout << " pitch: " <<  pitch << std::endl;
	  const double length = layergeom->get_pixel_z();
	  //	std::cout << " length: " << length << std::endl;
	  const double phisize = phibins.size() * pitch;
	  const double zsize = zbins.size() * length;
	  
	  static const double invsqrt12 = 1./std::sqrt(12);
	  
	  // scale factors (phi direction)
	  /*
	    they corresponds to clusters of size (2,2), (2,3), (3,2) and (3,3) in phi and z
	    other clusters, which are very few and pathological, get a scale factor of 1
	    These scale factors are applied to produce cluster pulls with width unity
	  */
	  
	  double phierror = pitch * invsqrt12;
	  
	  static constexpr std::array<double, 7> scalefactors_phi = {{ 0.36, 0.6,0.37,0.49,0.4,0.37,0.33 }};
	  if(phibins.size() == 1 && zbins.size() == 1) phierror*=scalefactors_phi[0];
	  else if(phibins.size() == 2 && zbins.size() == 1) phierror*=scalefactors_phi[1];
	  else if(phibins.size() == 1 && zbins.size() == 2) phierror*=scalefactors_phi[2];
	  else if( phibins.size() == 2 && zbins.size() == 2 ) phierror*=scalefactors_phi[0];
	  else if( phibins.size() == 2 && zbins.size() == 3 )  phierror*=scalefactors_phi[1];
	  else if( phibins.size() == 3 && zbins.size() == 2 )  phierror*=scalefactors_phi[2];
	  else if( phibins.size() == 3 && zbins.size() == 3 )  phierror*=scalefactors_phi[3];
	  
	  
	  // scale factors (z direction)
	  /*
	    they corresponds to clusters of size (2,2), (2,3), (3,2) and (3,3) in z and phi
	    other clusters, which are very few and pathological, get a scale factor of 1
	  */
	  static constexpr std::array<double, 4> scalefactors_z = {{ 0.47, 0.48, 0.71, 0.55 }};
	  double zerror = length*invsqrt12;
	  if( zbins.size() == 2 && phibins.size() == 2 ) zerror*=scalefactors_z[0];
	  else if( zbins.size() == 2 && phibins.size() == 3 )  zerror*=scalefactors_z[1];
	  else if( zbins.size() == 3 && phibins.size() == 2 )  zerror*=scalefactors_z[2];
	  else if( zbins.size() == 3 && phibins.size() == 3 )  zerror*=scalefactors_z[3];
	  
	  if(Verbosity() > 0)
	    cout << " MvtxClusterizer: cluskey " << ckey << " layer " << layer << " rad " << layergeom->get_radius() << " phibins " << phibins.size() << " pitch " << pitch << " phisize " << phisize 
		 << " zbins " << zbins.size() << " length " << length << " zsize " << zsize 
		 << " local x " << locclusx << " local y " << locclusz
		 << endl;
	  
	  if(m_cluster_version==3){
	    auto clus = std::make_unique<TrkrClusterv3>();
	    clus->setAdc(nhits);
	    clus->setLocalX(locclusx);
	    clus->setLocalY(locclusz);
	    // Take the rphi and z uncertainty of the cluster
	    clus->setActsLocalError(0,0,square(phierror));
	    clus->setActsLocalError(0,1,0.);
	    clus->setActsLocalError(1,0,0.);
	    clus->setActsLocalError(1,1,square(zerror));
	    
	    // All silicon surfaces have a 1-1 map to hitsetkey. 
	    // So set subsurface key to 0
	    clus->setSubSurfKey(0);
	    
	    if (Verbosity() > 2)
	      clus->identify();
	    
	    m_clusterlist->addClusterSpecifyKey(ckey, clus.release());
	  }else if(m_cluster_version==4){
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
	}  // clusitr loop
    }  // loop over hitsets
      
  if(Verbosity() > 1)
    {
      // check that the associations were written correctly
      m_clusterhitassoc->identify();
    }
      
  return;
}

void MvtxClusterizer::PrintClusters(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    TrkrClusterContainer *clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    if (!clusterlist) return;

    cout << "================= After MvtxClusterizer::process_event() ====================" << endl;

    cout << " There are " << clusterlist->size() << " clusters recorded: " << endl;

    if(Verbosity() > 3)  clusterlist->identify();

    cout << "===========================================================================" << endl;
  }

  return;
}

