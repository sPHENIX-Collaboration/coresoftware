#include "MvtxClusterizer.h"

#include <g4hough/SvtxHitMap.h>
#include <g4hough/SvtxHit.h>
#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxClusterMap_v1.h>
#include <g4hough/SvtxCluster.h>
#include <g4hough/SvtxCluster_v1.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeom_MAPS.h>
#include <g4detectors/PHG4CylinderGeom_Siladders.h>
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CylinderCellGeom.h>

#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>

#include <TMatrixF.h>
#include <TVector3.h>

#define BOOST_NO_HASH // Our version of boost.graph is incompatible with GCC-4.3 w/o this flag
#include <boost/bind.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
using namespace boost;

#include <iostream>
#include <stdexcept>
#include <cmath>

using namespace std;

static const float twopi = 2.0*M_PI;

bool MvtxClusterizer::maps_ladder_lessthan(const PHG4Cell* lhs, 
					  const PHG4Cell* rhs) 
{


    if( lhs->get_phibin() < rhs->get_phibin() ) return true;
    else if( lhs->get_phibin() == rhs->get_phibin() ){
      if( lhs->get_zbin() < rhs->get_zbin() ) return true;
    }
    
    
  return false;
}

bool MvtxClusterizer::maps_ladder_are_adjacent(const PHG4Cell* lhs, 
					      const PHG4Cell* rhs) 
{
  int lhs_layer = lhs->get_layer();
  int rhs_layer = rhs->get_layer();
  if (lhs_layer != rhs_layer) return false;

  // want to cluster only within a chip
  if(lhs->get_stave_index() != rhs->get_stave_index()) return false;
  if(lhs->get_half_stave_index() != rhs->get_half_stave_index()) return false;
  if(lhs->get_module_index() != rhs->get_module_index()) return false;
  if(lhs->get_chip_index() != rhs->get_chip_index()) return false;
  
  if (get_z_clustering(lhs_layer)) {
    if( fabs(lhs->get_zbin() - rhs->get_zbin()) <= 1 ) {
      if( fabs(lhs->get_phibin() - rhs->get_phibin()) <= 1 ){
	return true;
      }
    }
  } else {
    if( fabs(lhs->get_zbin() - rhs->get_zbin()) == 0 ) {
      if( fabs(lhs->get_phibin() - rhs->get_phibin()) <= 1 ){
	return true;
      } 
    }
  }

  return false;
}

MvtxClusterizer::MvtxClusterizer(const string &name,
					 unsigned int min_layer,
					 unsigned int max_layer) :
  SubsysReco(name),
  _hits(NULL),
  _clusterlist(NULL),
  _fraction_of_mip(0.5),
  _thresholds_by_layer(),
  _make_z_clustering(),
  _make_e_weights(),
  _min_layer(min_layer),
  _max_layer(max_layer),
  _timer(PHTimeServer::get()->insert_new(name)) 
{

}

int MvtxClusterizer::InitRun(PHCompositeNode* topNode) 
{

  // get node containing the digitized hits
  _hits = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
  if (!_hits) {
    cout << PHWHERE << "ERROR: Can't find node SvtxHitMap" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  //-----------------
  // Add Cluster Node
  //-----------------

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode 
    = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST"));
  if (!dstNode) {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter_dst(dstNode);
    
  // Create the SVX node if required
  PHCompositeNode* svxNode 
    = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode","SVTX"));
  if (!svxNode) {
    svxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svxNode);
  }
  
  // Create the Cluster node if required
  SvtxClusterMap *svxclusters 
    = findNode::getClass<SvtxClusterMap>(dstNode,"SvtxClusterMap");
  if (!svxclusters) {
    svxclusters = new SvtxClusterMap_v1();
    PHIODataNode<PHObject> *SvtxClusterMapNode =
      new PHIODataNode<PHObject>(svxclusters, "SvtxClusterMap", "PHObject");
    svxNode->addNode(SvtxClusterMapNode);
  }

  //----------------
  // Report Settings
  //----------------

  if (verbosity > 0) {
    cout << "====================== MvtxClusterizer::InitRun() =====================" << endl;
    cout << " Fraction of expected thickness MIP energy = " << _fraction_of_mip << endl;
    for (std::map<int,float>::iterator iter = _thresholds_by_layer.begin();
	 iter != _thresholds_by_layer.end();
	 ++iter) {
      cout << " Cluster Threshold in Layer #" << iter->first << " = " << 1.0e6*iter->second << " keV" << endl;
    }
    for (std::map<int,bool>::iterator iter = _make_z_clustering.begin();
	 iter != _make_z_clustering.end();
	 ++iter) {
      cout << " Z-dimension Clustering in Layer #" << iter->first << " = " << boolalpha << iter->second << noboolalpha << endl;
    }
    for (std::map<int,bool>::iterator iter = _make_e_weights.begin();
	 iter != _make_e_weights.end();
	 ++iter) {
      cout << " Energy weighting clusters in Layer #" << iter->first << " = " << boolalpha << iter->second << noboolalpha << endl;
    }
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxClusterizer::process_event(PHCompositeNode *topNode) {

  _timer.get()->restart();
  
  _clusterlist = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
  if (!_clusterlist) 
    {
      cout << PHWHERE << " ERROR: Can't find SvtxClusterMap." << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  _clusterlist->Reset();
  
  ClusterMapsLadderCells(topNode);

  PrintClusters(topNode);
  
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

void MvtxClusterizer::ClusterMapsLadderCells(PHCompositeNode *topNode) 
{

  if(verbosity > 0)
    cout << "Entering MvtxClusterizer::ClusterMapsLadderCells " << endl;

  //----------
  // Get Nodes
  //----------

  // get the SVX geometry object
  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode,"CYLINDERGEOM_MAPS");
  if (!geom_container) return;
  
  PHG4HitContainer* g4hits =  findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_MAPS");
  if (!g4hits) return;
  
  PHG4CellContainer* cells =  findNode::getClass<PHG4CellContainer>(topNode,"G4CELL_MAPS");
  if (!cells) return; 
 
  //-----------
  // Clustering
  //-----------

  // sort hits layer by layer
  std::multimap<int,SvtxHit*> layer_hits_mmap;
  for (SvtxHitMap::Iter iter = _hits->begin();
       iter != _hits->end();
       ++iter) {
    SvtxHit* hit = iter->second;
    layer_hits_mmap.insert(make_pair(hit->get_layer(),hit));
  }


  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) {

    int layer = layeriter->second->get_layer();

    if ((unsigned int)layer < _min_layer) continue;
    if ((unsigned int)layer > _max_layer) continue;

    std::map<PHG4Cell*,SvtxHit*> cell_hit_map;
    vector<PHG4Cell*> cell_list;
    for (std::multimap<int,SvtxHit*>::iterator hiter = layer_hits_mmap.lower_bound(layer);
	 hiter != layer_hits_mmap.upper_bound(layer);
	 ++hiter) {
      SvtxHit* hit = hiter->second;
      PHG4Cell* cell = cells->findCell(hit->get_cellid());
      cell_list.push_back(cell);
      cell_hit_map.insert(make_pair(cell,hit));
    }
    
    if (cell_list.size() == 0) continue; // if no cells, go to the next layer

    double hitx=0, hity=0, hitz =0;

    if(verbosity >4)
      {
	cout << "get g4_hit hit positions for layer  " << layer << endl;
	PHG4HitContainer::ConstIterator hiter;
	PHG4HitContainer::ConstRange hit_begin_end = g4hits->getHits(layer);
	for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
	  {    
	    hitx = hiter->second->get_x(0);
	    hity = hiter->second->get_y(0);
	    hitz = hiter->second->get_z(0);
	    cout << "  hit " << hiter->second->get_hit_id() << "  hitx " << hitx << " hity " << hity << " hitz " << hitz << endl;
	  }
      }
    
    // i'm not sure this sorting is ever really used
    sort(cell_list.begin(), cell_list.end(), MvtxClusterizer::maps_ladder_lessthan);

    typedef adjacency_list <vecS, vecS, undirectedS> Graph;
    Graph G;

    for(unsigned int i=0; i<cell_list.size(); i++) {
      for(unsigned int j=i+1; j<cell_list.size(); j++) {
        if(maps_ladder_are_adjacent(cell_list[i], cell_list[j]) )
          add_edge(i,j,G);
      }
      
      add_edge(i,i,G);
    }


    // Find the connections between the vertices of the graph (vertices are the rawhits, 
    // connections are made when they are adjacent to one another)
    vector<int> component(num_vertices(G));
    
    // this is the actual clustering, performed by boost
    connected_components(G, &component[0]); 

    // Loop over the components(hit cells) compiling a list of the
    // unique connected groups (ie. clusters).
    set<int> cluster_ids; // unique components
    multimap<int, PHG4Cell*> clusters;
    for (unsigned int i=0; i<component.size(); i++) {
      cluster_ids.insert( component[i] );
      clusters.insert( make_pair(component[i], cell_list[i]) );
    }
    
    // 
    for (set<int>::iterator clusiter = cluster_ids.begin(); 
	 clusiter != cluster_ids.end(); 
	 clusiter++ ) {
      
      int clusid = *clusiter;
      pair<multimap<int, PHG4Cell*>::iterator,
	   multimap<int, PHG4Cell*>::iterator> clusrange = clusters.equal_range(clusid);
      
      multimap<int, PHG4Cell*>::iterator mapiter = clusrange.first;
      
      int layer = mapiter->second->get_layer();
      PHG4CylinderGeom_MAPS *geom = (PHG4CylinderGeom_MAPS*) geom_container->GetLayerGeom(layer);

      if(verbosity > 2)
	cout << "Filling cluster id " << clusid << " in  layer " << layer << endl;
      
      SvtxCluster_v1 clus;
      clus.set_layer( layer );
      float clus_energy = 0.0;
      unsigned int clus_adc = 0;

      // determine the size of the cluster in phi and z
      // useful for track fitting the cluster

      set<int> phibins;
      set<int> zbins;
      for (mapiter = clusrange.first; mapiter != clusrange.second; mapiter++ ) {
	PHG4Cell* cell = mapiter->second;     

	int pixel_number = cell->get_pixel_index();
	// binphi is the cell index in the phi direction in the sensor
	int binphi = geom->get_pixel_X_from_pixel_number(pixel_number);
	phibins.insert(binphi);
	// binz is the cell index in the z direction in the sensor
	int binz = geom->get_pixel_Z_from_pixel_number(pixel_number);
	zbins.insert(binz);

	if(verbosity > 2)
	  cout << "   pixel number " << pixel_number << " binphi = " << binphi  << " binz = " << binz  << endl;
      }

      float thickness = geom->get_pixel_thickness();
      float pitch = geom->get_pixel_x();
      float length = geom->get_pixel_z();
      float phisize = phibins.size()*pitch;
      float zsize = zbins.size()*length;
      // tilt refers to a rotation around the radial vector from the origin, and this is zero for the MAPS ladders
      float tilt = 0.0;

      // determine the cluster position...
      double xsum = 0.0;
      double ysum = 0.0;
      double zsum = 0.0;
      unsigned nhits = 0;

      int stave_index = -1;
      int half_stave_index = -1;
      int module_index = -1;
      int chip_index = -1;

      for(mapiter = clusrange.first; mapiter != clusrange.second; mapiter++ ) {
        PHG4Cell* cell = mapiter->second;

	//if(verbosity > 2)	
	  //cell->identify();
	
	SvtxHit* hit = cell_hit_map[cell];
	
	clus.insert_hit(hit->get_id());
	
        clus_energy += hit->get_e();
	clus_adc    += hit->get_adc();

	// find the center of the pixel in world coords
	TVector3 local_coords = geom->get_local_coords_from_pixel(cell->get_pixel_index());
	TVector3 world_coords = geom->get_world_from_local_coords(cell->get_stave_index(), cell->get_half_stave_index(), cell->get_module_index(), cell->get_chip_index(), local_coords);
	double hit_location[3] = {world_coords.X(), world_coords.Y(), world_coords.Z()};

	// These will be used later to get the sensor position so that the sensor phi can be calculated
	stave_index = cell->get_stave_index();
	half_stave_index = cell->get_half_stave_index();
	module_index = cell->get_module_index();
	chip_index = cell->get_chip_index();

	if (_make_e_weights[layer]) {
	  xsum += hit_location[0] * hit->get_adc();
	  ysum += hit_location[1] * hit->get_adc();
	  zsum += hit_location[2] * hit->get_adc();  
	} else {
	  xsum += hit_location[0];
	  ysum += hit_location[1];
	  zsum += hit_location[2];
	}
	
	if(verbosity > 2)
	  {
	    cout << "  From  geometry object: hit x " << hit_location[0] << " hit y " << hit_location[1] << " hit z " << hit_location[2] << " from hit object: e " << hit->get_e() << " hit adc " << hit->get_adc() << " e weight " << _make_e_weights[layer] << endl;
	  }
	
	++nhits;
      }
      
      double clusx = NAN;
      double clusy = NAN;
      double clusz = NAN;

      if (_make_e_weights[layer]) {
	clusx = xsum / clus_adc;
	clusy = ysum / clus_adc;
	clusz = zsum / clus_adc;	
      } else {
	clusx = xsum / nhits;
	clusy = ysum / nhits;
	clusz = zsum / nhits;
      }
      
      double ladder_location[3] = {0.0,0.0,0.0};
      // returns the center of the sensor in world coordinates - used to get the ladder phi location
      geom->find_sensor_center(stave_index, half_stave_index, module_index, chip_index, ladder_location);
      double ladderphi = atan2( ladder_location[1], ladder_location[0] );
      ladderphi += geom->get_stave_phi_tilt();

      //cout << "sensor center = " << ladder_location[0] << " " << ladder_location[1] << " " << ladder_location[2] << endl;            

      clus.set_position(0, clusx);
      clus.set_position(1, clusy);
      clus.set_position(2, clusz);

      clus.set_e(clus_energy);
      clus.set_adc(clus_adc);

      float invsqrt12 = 1.0/sqrt(12.0);
      
      TMatrixF DIM(3,3);
      DIM[0][0] = pow(0.5*thickness,2);
      DIM[0][1] = 0.0;
      DIM[0][2] = 0.0;
      DIM[1][0] = 0.0;
      DIM[1][1] = pow(0.5*phisize,2);
      DIM[1][2] = 0.0;
      DIM[2][0] = 0.0;
      DIM[2][1] = 0.0;
      DIM[2][2] = pow(0.5*zsize,2);

      TMatrixF ERR(3,3);
      ERR[0][0] = pow(0.5*thickness*invsqrt12,2);
      ERR[0][1] = 0.0;
      ERR[0][2] = 0.0;
      ERR[1][0] = 0.0;
      ERR[1][1] = pow(0.5*phisize*invsqrt12,2);
      ERR[1][2] = 0.0;
      ERR[2][0] = 0.0;
      ERR[2][1] = 0.0;
      ERR[2][2] = pow(0.5*zsize*invsqrt12,2);
      
      TMatrixF ROT(3,3);
      ROT[0][0] = cos(ladderphi);
      ROT[0][1] = -1.0*sin(ladderphi);
      ROT[0][2] = 0.0;
      ROT[1][0] = sin(ladderphi);
      ROT[1][1] = cos(ladderphi);
      ROT[1][2] = 0.0;
      ROT[2][0] = 0.0;
      ROT[2][1] = 0.0;
      ROT[2][2] = 1.0;

      TMatrixF TILT(3,3);
      TILT[0][0] = 1.0;
      TILT[0][1] = 0.0;
      TILT[0][2] = 0.0;
      TILT[1][0] = 0.0;
      TILT[1][1] = cos(tilt);
      TILT[1][2] = -1.0*sin(tilt);
      TILT[2][0] = 0.0;
      TILT[2][1] = sin(tilt);
      TILT[2][2] = cos(tilt);

      TMatrixF R(3,3);
      R = ROT * TILT;
      
      TMatrixF R_T(3,3);
      R_T.Transpose(R);
          
      TMatrixF COVAR_DIM(3,3);
      COVAR_DIM = R * DIM * R_T;
      
      clus.set_size( 0 , 0 , COVAR_DIM[0][0] );
      clus.set_size( 0 , 1 , COVAR_DIM[0][1] );
      clus.set_size( 0 , 2 , COVAR_DIM[0][2] );
      clus.set_size( 1 , 0 , COVAR_DIM[1][0] );
      clus.set_size( 1 , 1 , COVAR_DIM[1][1] );
      clus.set_size( 1 , 2 , COVAR_DIM[1][2] );
      clus.set_size( 2 , 0 , COVAR_DIM[2][0] );
      clus.set_size( 2 , 1 , COVAR_DIM[2][1] );
      clus.set_size( 2 , 2 , COVAR_DIM[2][2] );

      TMatrixF COVAR_ERR(3,3);
      COVAR_ERR = R * ERR * R_T;
      
      clus.set_error( 0 , 0 , COVAR_ERR[0][0] );
      clus.set_error( 0 , 1 , COVAR_ERR[0][1] );
      clus.set_error( 0 , 2 , COVAR_ERR[0][2] );
      clus.set_error( 1 , 0 , COVAR_ERR[1][0] );
      clus.set_error( 1 , 1 , COVAR_ERR[1][1] );
      clus.set_error( 1 , 2 , COVAR_ERR[1][2] );
      clus.set_error( 2 , 0 , COVAR_ERR[2][0] );
      clus.set_error( 2 , 1 , COVAR_ERR[2][1] );
      clus.set_error( 2 , 2 , COVAR_ERR[2][2] );
      
      // if (clus_energy > get_threshold_by_layer(layer)) {
      if (true) {
	SvtxCluster* ptr = _clusterlist->insert(&clus);
	if (!ptr->isValid()) {
	  static bool first = true;
	  if (first) {
	    cout << PHWHERE << "ERROR: Invalid SvtxClusters are being produced" << endl;
	    ptr->identify();
	    first = false;
	  }
	}
	if(verbosity > 1)
	  {
	    // fairly complete sanity check:
	    // Get the list of g4hit positions for this cluster and compare positions	    
	    cout << " For cluster " << ptr->get_id() << endl;
	    cout << " cluster position: x " << ptr->get_x() << " y " << ptr->get_y() << " z " << ptr->get_z() << endl;
	    cout << " list of hit id's: " << endl;
	    for (SvtxCluster::HitIter iter = ptr->begin_hits(); iter != ptr->end_hits(); ++iter) {
	      cout << "  " << *iter << " ";
	      SvtxHit *hit = _hits->get(*iter);
	      cout << " cell id from hit = " << hit->get_cellid() << " : " << endl; 
	      PHG4Cell *cell = cells->findCell(hit->get_cellid());
	      cout  << " cell data: stave " << cell->get_property_int(PHG4Cell::prop_stave_index)
		    << " half stave " << cell->get_property_int(PHG4Cell::prop_half_stave_index)
		    << " module " << cell->get_property_int(PHG4Cell::prop_module_index)
		    << " chip " << cell->get_property_int(PHG4Cell::prop_chip_index)
		    << " pixel_index " << cell->get_pixel_index()
		    << " phibin " << cell->get_phibin()
		    << " zbin " << cell->get_zbin()
		    << endl;
	      
	      for (PHG4Cell::EdepConstIterator g4iter = cell->get_g4hits().first;
		   g4iter != cell->get_g4hits().second;
		   ++g4iter) {		
		PHG4Hit *g4hit = g4hits->findHit(g4iter->first);
		cout << "    g4hit position: x " << g4hit->get_x(0) << " y " << g4hit->get_y(0) << " z " << g4hit->get_z(0) << endl;		

		// test that there is not a large difference between the cluster position and the position(s) of the g4 hits that contributed to it
		if( fabs(ptr->get_z() - g4hit->get_z(0)) > 0.1)
		  cout << "       ALERT! g4hit entry point and cluster Z do not agree by " << endl;
	      }
	    }
	    cout << endl;
	  }

	if (verbosity>1) {
	  double radius = sqrt(clusx*clusx+clusy*clusy);
	  double clusphi = atan2(clusy,clusx);
	  cout << "clus energy " << clus_energy << " clus_adc " << clus_adc << " r=" << radius << " phi=" << clusphi << " z=" << clusz << endl;
	  cout << "pos=(" << clus.get_position(0) << ", " << clus.get_position(1)
	       << ", " << clus.get_position(2) << ")" << endl;
	  cout << endl;
	}
      }	else if (verbosity>1) {
	double radius = sqrt(clusx*clusx+clusy*clusy);
	double clusphi = atan2(clusy,clusx);
	cout << "MAPS ladder cell: removed, clus_energy = " << clus_energy << " below threshold of " <<  0  << " clus_adc " << clus_adc <<  " r=" << radius << " phi=" << clusphi << " z=" << clusz << endl;
	cout << "pos=(" << clus.get_position(0) << ", " << clus.get_position(1)
	     << ", " << clus.get_position(2) << ")" << endl;
	cout << endl;
      } 
    }
  }
  
  return;
}


void MvtxClusterizer::PrintClusters(PHCompositeNode *topNode) {

  if (verbosity >= 1) {

    SvtxClusterMap *clusterlist = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
    if (!clusterlist) return;
    
    cout << "================= MvtxClusterizer::process_event() ====================" << endl;
  

    cout << " Found and recorded the following " << clusterlist->size() << " clusters: " << endl;

    unsigned int icluster = 0;
    for (SvtxClusterMap::Iter iter = clusterlist->begin();
	 iter != clusterlist->end();
	 ++iter) {

      SvtxCluster* cluster = iter->second;
      cout << icluster << " of " << clusterlist->size() << endl;
      cluster->identify();
      ++icluster;
    }
    
    cout << "===========================================================================" << endl;
  }
  
  return;
}
