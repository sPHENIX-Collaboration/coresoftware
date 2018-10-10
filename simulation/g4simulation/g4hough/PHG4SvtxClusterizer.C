#include "PHG4SvtxClusterizer.h"

#include "SvtxHitMap.h"
#include "SvtxHit.h"
#include "SvtxClusterMap.h"
#include "SvtxClusterMap_v1.h"
#include "SvtxCluster.h"
#include "SvtxCluster_v1.h"

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
#include <g4detectors/PHG4CylinderGeomSiLadders.h>
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

bool PHG4SvtxClusterizer::lessthan(const PHG4Cell* lhs, 
				   const PHG4Cell* rhs) {
  int lhsphibin = PHG4CellDefs::SizeBinning::get_phibin(lhs->get_cellid());
  int rhsphibin =  PHG4CellDefs::SizeBinning::get_phibin(rhs->get_cellid());
  int lhszbin = PHG4CellDefs::SizeBinning::get_zbin(lhs->get_cellid());
  int rhszbin = PHG4CellDefs::SizeBinning::get_zbin(rhs->get_cellid());

  if( lhsphibin < rhsphibin ) return true;
  else if( lhsphibin == rhsphibin ){
    if( lhszbin < rhszbin ) return true;
  }

  return false;
}

bool PHG4SvtxClusterizer::ladder_lessthan(const PHG4Cell* lhs, 
					  const PHG4Cell* rhs) {

  if ( lhs->get_ladder_z_index() == rhs->get_ladder_z_index() &&
       lhs->get_ladder_phi_index() == rhs->get_ladder_phi_index())
 { 

    if( lhs->get_phibin() < rhs->get_phibin() ) return true;
    else if( lhs->get_phibin() == rhs->get_phibin() ){
      if( lhs->get_zbin() < rhs->get_zbin() ) return true;
    }
    
  } else {
    if ( lhs->get_zbin() < rhs->get_zbin() ) return true;   
  }
    
  return false;
}

bool PHG4SvtxClusterizer::maps_ladder_lessthan(const PHG4Cell* lhs, 
					  const PHG4Cell* rhs) {


    if( lhs->get_phibin() < rhs->get_phibin() ) return true;
    else if( lhs->get_phibin() == rhs->get_phibin() ){
      if( lhs->get_zbin() < rhs->get_zbin() ) return true;
    }
    
    
  return false;
}

bool PHG4SvtxClusterizer::are_adjacent(const PHG4Cell* lhs, 
				       const PHG4Cell* rhs,
                                       const int &nphibins) {

  int lhs_layer = lhs->get_layer();
  int rhs_layer = rhs->get_layer();
  if (lhs_layer != rhs_layer) return false;

  int lhsphibin = PHG4CellDefs::SizeBinning::get_phibin(lhs->get_cellid());
  int rhsphibin =  PHG4CellDefs::SizeBinning::get_phibin(rhs->get_cellid());
  int lhszbin = PHG4CellDefs::SizeBinning::get_zbin(lhs->get_cellid());
  int rhszbin = PHG4CellDefs::SizeBinning::get_zbin(rhs->get_cellid());
  if (get_z_clustering(lhs_layer)) {
    if( fabs(lhszbin - rhszbin) <= 1 ) {
      if( fabs(lhsphibin - rhsphibin) <= 1 ){
	return true;
      } else if(lhsphibin == 0 || rhsphibin == 0) {
	if(fabs(lhsphibin - rhsphibin) == (nphibins-1))
	  return true;
      }
    }
  } else {
    if( fabs(lhszbin - rhszbin) == 0 ) {
      if( fabs(lhsphibin - rhsphibin) <= 1 ){
	return true;
      } else if(lhsphibin == 0 || rhsphibin == 0) {
	if(fabs(lhsphibin - rhsphibin) == (nphibins-1))
	  return true;
      }
    }
  }

  return false;
}

bool PHG4SvtxClusterizer::maps_ladder_are_adjacent(const PHG4Cell* lhs, 
					      const PHG4Cell* rhs) {
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

bool PHG4SvtxClusterizer::ladder_are_adjacent(const PHG4Cell* lhs, const PHG4Cell* rhs) {

  if(Verbosity() > 2)
    {
      cout << " lhs layer " <<  lhs->get_layer() 
	   << " lhs ladder z index " <<  lhs->get_ladder_z_index() 
	   << " lhs ladder phi index " <<  lhs->get_ladder_phi_index()
	   << " lhs z bin " <<  lhs->get_zbin()
	   << " lhs phi bin " <<  lhs->get_phibin()
	   << endl;
      
      cout << " rhs layer " <<  rhs->get_layer() 
	   << " rhs ladder z index " <<  rhs->get_ladder_z_index() 
	   << " rhs ladder phi index " <<  rhs->get_ladder_phi_index()
	   << " rhs z bin " <<  rhs->get_zbin()
	   << " rhs phi bin " <<  rhs->get_phibin()
	   << endl;
    }

  int lhs_layer = lhs->get_layer();
  int rhs_layer = rhs->get_layer();
  if (lhs_layer != rhs_layer) return false;

  if ( !( lhs->get_ladder_z_index() == rhs->get_ladder_z_index() &&
	  lhs->get_ladder_phi_index() == rhs->get_ladder_phi_index()) ) return false;
  
  if (get_z_clustering(lhs_layer)) {
    if( fabs(lhs->get_zbin() - rhs->get_zbin()) <= 1 ) {
      if( fabs(lhs->get_phibin() - rhs->get_phibin()) <= 1 ){
	return true;
      }
    }
  } else {
    if( fabs(lhs->get_zbin() - rhs->get_zbin()) == 0 ) {
      if( fabs(lhs->get_phibin() - rhs->get_phibin()) <= 1 ){
	if(Verbosity() > 2) cout << "    accepted " << endl;
	return true;
      } 
    }
  }

  return false;
}

PHG4SvtxClusterizer::PHG4SvtxClusterizer(const string &name,
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
  _timer(PHTimeServer::get()->insert_new(name)) {}

int PHG4SvtxClusterizer::InitRun(PHCompositeNode* topNode) {

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

  //---------------------
  // Calculate Thresholds
  //---------------------
  
  CalculateCylinderThresholds(topNode);
  CalculateLadderThresholds(topNode);
  CalculateMapsLadderThresholds(topNode);

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0) {
    cout << "====================== PHG4SvtxClusterizer::InitRun() =====================" << endl;
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

int PHG4SvtxClusterizer::process_event(PHCompositeNode *topNode) {

  _timer.get()->restart();
  
  _clusterlist = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
  if (!_clusterlist) 
    {
      cout << PHWHERE << " ERROR: Can't find SvtxClusterMap." << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  _clusterlist->Reset();
  
  ClusterCylinderCells(topNode);
  ClusterLadderCells(topNode);
  ClusterMapsLadderCells(topNode);

  PrintClusters(topNode);
  
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4SvtxClusterizer::CalculateCylinderThresholds(PHCompositeNode *topNode) {

  // get the SVX geometry object
  PHG4CylinderCellGeomContainer* geom_container = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode,"CYLINDERCELLGEOM_SVTX");
  if (!geom_container) return;
  
  // determine cluster thresholds and layer index mapping
  PHG4CylinderCellGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) {

    // index mapping
    int layer = layeriter->second->get_layer();

    // cluster threshold
    float thickness = (layeriter->second)->get_thickness();
    float threshold = _fraction_of_mip*0.003876*thickness;
    _thresholds_by_layer.insert(std::make_pair(layer,threshold));

    // fill in a default z_clustering value if not present
    if (_make_z_clustering.find(layer) == _make_z_clustering.end()) {
      _make_z_clustering.insert(std::make_pair(layer,true));
    }
    
    if (_make_e_weights.find(layer) == _make_e_weights.end()) {
      _make_e_weights.insert(std::make_pair(layer,false));
    }
  }
  
  return;
}

void PHG4SvtxClusterizer::CalculateLadderThresholds(PHCompositeNode *topNode) {

  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode,"G4CELL_SILICON_TRACKER");
  if (!cells) return;

  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode,"CYLINDERGEOM_SILICON_TRACKER");
  if (!geom_container) return;
  
  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) {

    // index mapping
    int layer = layeriter->second->get_layer();

    // cluster threshold
    float thickness = (layeriter->second)->get_thickness();
    float threshold = _fraction_of_mip*0.003876*thickness;
    _thresholds_by_layer.insert(std::make_pair(layer,threshold));

    // fill in a default z_clustering value if not present
    if (_make_z_clustering.find(layer) == _make_z_clustering.end()) {
      _make_z_clustering.insert(std::make_pair(layer,true));
    }

    if (_make_e_weights.find(layer) == _make_e_weights.end()) {
      _make_e_weights.insert(std::make_pair(layer,false));
    }
  }
  
  return;
}

void PHG4SvtxClusterizer::CalculateMapsLadderThresholds(PHCompositeNode *topNode) {

  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode,"G4CELL_MAPS");
  if (!cells) return;

  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode,"CYLINDERGEOM_MAPS");
  if (!geom_container) return;
  
  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) {

    // index mapping
    int layer = layeriter->second->get_layer();

    // cluster threshold
    float thickness = (layeriter->second)->get_pixel_thickness();
    float threshold = _fraction_of_mip*0.003876*thickness;
    _thresholds_by_layer.insert(std::make_pair(layer,threshold));

    // fill in a default z_clustering value if not present
    if (_make_z_clustering.find(layer) == _make_z_clustering.end()) {
      _make_z_clustering.insert(std::make_pair(layer,true));
    }

    if (_make_e_weights.find(layer) == _make_e_weights.end()) {
      _make_e_weights.insert(std::make_pair(layer,false));
    }
  }
  
  return;
}

void PHG4SvtxClusterizer::ClusterCylinderCells(PHCompositeNode *topNode) {

  //----------
  // Get Nodes
  //----------

  // get the SVX geometry object
  PHG4CylinderCellGeomContainer* geom_container = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode,"CYLINDERCELLGEOM_SVTX");
  if (!geom_container) return;
  
  PHG4HitContainer* g4hits = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_SVTX");
  if (!g4hits) return;
  
  PHG4CellContainer* cells = findNode::getClass<PHG4CellContainer>(topNode,"G4CELL_SVTX");
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
  
  // loop over cylinder layers
  PHG4CylinderCellGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) {

    int layer = layeriter->second->get_layer();

    if ((unsigned int)layer < _min_layer) continue;
    if ((unsigned int)layer > _max_layer) continue;
    
    int nphibins = layeriter->second->get_phibins();

    // loop over all hits/cells in this layer
    std::map<PHG4Cell*,SvtxHit*> cell_hit_map;
    std::vector<PHG4Cell*> cell_list;   
    for (std::multimap<int,SvtxHit*>::iterator hiter = layer_hits_mmap.lower_bound(layer);
	 hiter != layer_hits_mmap.upper_bound(layer);
	 ++hiter) {
      SvtxHit* hit = hiter->second;
      PHG4Cell* cell = cells->findCell(hit->get_cellid());
      cell_list.push_back(cell);
      cell_hit_map.insert(make_pair(cell,hit));
    }

    if (cell_list.size() == 0) continue; // if no cells, go to the next layer
    
    sort(cell_list.begin(), cell_list.end(), PHG4SvtxClusterizer::lessthan);

    typedef adjacency_list <vecS, vecS, undirectedS> Graph;
    Graph G;

    for(unsigned int i=0; i<cell_list.size(); i++) {
      for(unsigned int j=i+1; j<cell_list.size(); j++) {
        if( are_adjacent(cell_list[i], cell_list[j], nphibins) )
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
    
    typedef multimap<int, PHG4Cell*>::iterator mapiterator;
    
    for (set<int>::iterator clusiter = cluster_ids.begin(); 
	 clusiter != cluster_ids.end(); 
	 clusiter++ ) {
      
      int clusid = *clusiter;
      pair<mapiterator,mapiterator> clusrange = clusters.equal_range(clusid);
      
      mapiterator mapiter = clusrange.first;
      
      int layer = mapiter->second->get_layer();
      PHG4CylinderCellGeom* geom = geom_container->GetLayerCellGeom(layer);
      
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
	
	phibins.insert(PHG4CellDefs::SizeBinning::get_phibin(cell->get_cellid()));
	zbins.insert(PHG4CellDefs::SizeBinning::get_zbin(cell->get_cellid()));
      }
      
      float pitch = geom->get_phistep()*(geom->get_radius()+0.5*geom->get_thickness());
      float thickness = geom->get_thickness();
      float length = geom->get_zstep();
      float phisize = phibins.size()*pitch;
      float zsize = zbins.size()*length;

      double xsum = 0.0;
      double ysum = 0.0;
      double zsum = 0.0;
      unsigned int nhits = 0;

      for(mapiter = clusrange.first; mapiter != clusrange.second; mapiter++ ) {
        PHG4Cell* cell = mapiter->second;
	SvtxHit* hit = cell_hit_map[cell];
	
	clus.insert_hit(hit->get_id());
	
        clus_energy += hit->get_e();
	clus_adc    += hit->get_adc();

	// compute the hit center
	double r   = geom->get_radius()+0.5*geom->get_thickness();
        double phi = geom->get_phicenter(PHG4CellDefs::SizeBinning::get_phibin(cell->get_cellid()));

	double x = r*cos(phi);
	double y = r*sin(phi);
        double z = geom->get_zcenter(PHG4CellDefs::SizeBinning::get_zbin(cell->get_cellid()));

	if (_make_e_weights[layer]) {
	  xsum += x * hit->get_adc();
	  ysum += y * hit->get_adc();
	  zsum += z * hit->get_adc();  
	} else {
	  xsum += x;
	  ysum += y;
	  zsum += z;
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
      
      double radius  = sqrt(clusx*clusx+clusy*clusy);
      double clusphi = atan2( clusy, clusx);
       
      clus.set_position( 0 , clusx );
      clus.set_position( 1 , clusy );
      clus.set_position( 2 , clusz );

      clus.set_e(clus_energy);
      clus.set_adc(clus_adc);

      float invsqrt12 = 1.0/sqrt(12.);
      
      TMatrixF DIM(3,3);
      DIM[0][0] = pow(0.0*0.5*thickness,2);
      DIM[0][1] = 0.0;
      DIM[0][2] = 0.0;
      DIM[1][0] = 0.0;
      DIM[1][1] = pow(0.5*phisize,2);
      DIM[1][2] = 0.0;
      DIM[2][0] = 0.0;
      DIM[2][1] = 0.0;
      DIM[2][2] = pow(0.5*zsize,2);

      const float corr_factor = 1.0; // cylinder

      TMatrixF ERR(3,3);
      ERR[0][0] = pow(0.0*thickness*invsqrt12*corr_factor,2);
      ERR[0][1] = 0.0;
      ERR[0][2] = 0.0;
      ERR[1][0] = 0.0;
      ERR[1][1] = pow(phisize*invsqrt12*corr_factor,2);
      ERR[1][2] = 0.0;
      ERR[2][0] = 0.0;
      ERR[2][1] = 0.0;
      ERR[2][2] = pow(zsize*invsqrt12*corr_factor,2);

      TMatrixF ROT(3,3);
      ROT[0][0] = cos(clusphi);
      ROT[0][1] = -sin(clusphi);
      ROT[0][2] = 0.0;
      ROT[1][0] = sin(clusphi);
      ROT[1][1] = cos(clusphi);
      ROT[1][2] = 0.0;
      ROT[2][0] = 0.0;
      ROT[2][1] = 0.0;
      ROT[2][2] = 1.0;

      TMatrixF ROT_T(3,3);
      ROT_T.Transpose(ROT);
      
      TMatrixF COVAR_DIM(3,3);
      COVAR_DIM = ROT * DIM * ROT_T;
      
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
      COVAR_ERR = ROT * ERR * ROT_T;

      clus.set_error( 0 , 0 , COVAR_ERR[0][0] );
      clus.set_error( 0 , 1 , COVAR_ERR[0][1] );
      clus.set_error( 0 , 2 , COVAR_ERR[0][2] );
      clus.set_error( 1 , 0 , COVAR_ERR[1][0] );
      clus.set_error( 1 , 1 , COVAR_ERR[1][1] );
      clus.set_error( 1 , 2 , COVAR_ERR[1][2] );
      clus.set_error( 2 , 0 , COVAR_ERR[2][0] );
      clus.set_error( 2 , 1 , COVAR_ERR[2][1] );
      clus.set_error( 2 , 2 , COVAR_ERR[2][2] );
      
      if (clus_energy > get_threshold_by_layer(layer)) {
	SvtxCluster* ptr = _clusterlist->insert(&clus);
	if (!ptr->isValid()) {
	  static bool first = true;
	  if (first) {
	    cout << PHWHERE << "ERROR: Invalid SvtxClusters are being produced" << endl;
	    ptr->identify();
	    first = false;
	  }
	}
	
	if (Verbosity()>1) {
	  cout << "r=" << radius << " phi=" << clusphi << " z=" << clusz << endl;
	  cout << "pos=(" << clus.get_position(0) << ", " << clus.get_position(1)
	       << ", " << clus.get_position(2) << ")" << endl;
	  cout << endl;
	}
      }	else if (Verbosity()>1) {
	cout << "silicon cylinder cell: removed, clus_energy = " << clus_energy << " below threshold of " <<  get_threshold_by_layer(layer)  << " clus_adc " << clus_adc <<  " r=" << radius << " phi=" << clusphi << " z=" << clusz << endl;
	cout << "pos=(" << clus.get_position(0) << ", " << clus.get_position(1)
	     << ", " << clus.get_position(2) << ")" << endl;
	cout << endl;
      } 
    }
  }
  
  return;
}

void PHG4SvtxClusterizer::ClusterLadderCells(PHCompositeNode *topNode) {

  if(Verbosity() > 0)
    cout << "Entering PHG4SvtxClusterizer::ClusterLadderCells " << endl;

  //----------
  // Get Nodes
  //----------

  // get the SVX geometry object
  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode,"CYLINDERGEOM_SILICON_TRACKER");
  if (!geom_container) return;
  
  PHG4HitContainer* g4hits =  findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_SILICON_TRACKER");
  if (!g4hits) return;
  
  PHG4CellContainer* cells =  findNode::getClass<PHG4CellContainer>(topNode,"G4CELL_SILICON_TRACKER");
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

  // this does nothing!
  for (std::multimap<int,SvtxHit*>::iterator it=layer_hits_mmap.begin(); it!=layer_hits_mmap.end(); ++it)
    {
      if( (*it).first < (int) _min_layer || (*it).first > (int) _max_layer) 
	continue;
    }

  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) {

    int layer = layeriter->second->get_layer();
    if(Verbosity() > 0)
      cout << " layer loop, current layer = " << layer << endl;

    if ((unsigned int)layer < _min_layer) continue;
    if ((unsigned int)layer > _max_layer) continue;

    std::map<PHG4Cell*,SvtxHit*> cell_hit_map;
    vector<PHG4Cell*> cell_list;
    for (std::multimap<int,SvtxHit*>::iterator hiter = layer_hits_mmap.lower_bound(layer);
	 hiter != layer_hits_mmap.upper_bound(layer);
	 ++hiter) {
      SvtxHit* hit = hiter->second;
      PHG4Cell* cell = cells->findCell(hit->get_cellid());
      if(Verbosity() > 2) 
	{
	  cout << "adding cell to cell_hit_map: ";
	  cell->print();
	}
      cell_list.push_back(cell);
      cell_hit_map.insert(make_pair(cell,hit));
    }

    if (cell_list.size() == 0) continue; // if no cells, go to the next layer
    
    // i'm not sure this sorting is ever really used
    sort(cell_list.begin(), cell_list.end(), PHG4SvtxClusterizer::ladder_lessthan);

    typedef adjacency_list <vecS, vecS, undirectedS> Graph;
    Graph G;

    // Find adjacent cell
    if(Verbosity() > 2) cout << "Find adjacent cells for layer " << layer << endl;
    for(unsigned int i=0; i<cell_list.size(); i++) {
      for(unsigned int j=i+1; j<cell_list.size(); j++) {
	if(Verbosity() > 2) 
	  {
	    cout << "compare cells " << i << " and " << j << endl;
	    cell_list[i]->print();
	    cell_list[j]->print();
	  }
        if(ladder_are_adjacent(cell_list[i], cell_list[j]) )
	  {
	    add_edge(i,j,G);
	    if(Verbosity() > 2) 
	      {
		cout << "Found edge " << i << "   " << j << endl;
		cell_list[i]->print();
		cell_list[j]->print();
	      }
	  }
      }
      
      add_edge(i,i,G);
    }
    if(Verbosity() > 2) cout << "finished looking for adjacent cells for layer " << layer << endl;

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
      //PHG4CylinderGeom* geom = geom_container->GetLayerGeom(layer);
      PHG4CylinderGeomSiLadders* geom = dynamic_cast<PHG4CylinderGeomSiLadders*> (geom_container->GetLayerGeom(layer));
      
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
	
	phibins.insert(cell->get_phibin());
	zbins.insert(cell->get_zbin());
      }

      float thickness = geom->get_thickness();
      float pitch = geom->get_strip_y_spacing();
      float length = geom->get_strip_z_spacing();
      float phisize = phibins.size()*pitch;
      float zsize = zbins.size()*length;

      // tilt refers to a rotation around the radial vector from the origin, and this is zero for the INTT ladders
      float tilt = 0; //geom->get_strip_tilt();

      // determine the cluster position...
      double xsum = 0.0;
      double ysum = 0.0;
      double zsum = 0.0;
      unsigned nhits = 0;

      int ladder_z_index = -1;
      int ladder_phi_index = -1;
      
      for(mapiter = clusrange.first; mapiter != clusrange.second; mapiter++ ) {
        PHG4Cell* cell = mapiter->second;
	SvtxHit* hit = cell_hit_map[cell];
	if(Verbosity()>0)
	  {
	    cout << "Add hit: "; 
	    hit->identify();
	    cout << " cell is ";
	    cell->print();
	  }
	clus.insert_hit(hit->get_id());
	
        clus_energy += hit->get_e();
	clus_adc    += hit->get_adc();

	double hit_location[3] = {0.0,0.0,0.0};
	geom->find_strip_center(cell->get_ladder_z_index(),
				cell->get_ladder_phi_index(),
				cell->get_zbin(),
				cell->get_phibin(),
				hit_location);

	ladder_z_index = cell->get_ladder_z_index();
	ladder_phi_index = cell->get_ladder_phi_index();

	if (_make_e_weights[layer]) {
	  xsum += hit_location[0] * hit->get_adc();
	  ysum += hit_location[1] * hit->get_adc();
	  zsum += hit_location[2] * hit->get_adc();  
	} else {
	  xsum += hit_location[0];
	  ysum += hit_location[1];
	  zsum += hit_location[2];
	}
	++nhits;
	if(Verbosity() > 2) cout << "     nhits = " << nhits << endl;
	if(Verbosity() > 2)
	  {
	    cout << "  From  geometry object: hit x " << hit_location[0] << " hit y " << hit_location[1] << " hit z " << hit_location[2]  << endl;
	    cout << "     nhits " << nhits << " clusx  = " << xsum/nhits << " clusy " << ysum/nhits << " clusz " << zsum/nhits  << endl;
	    
	    if( fabs(xsum/nhits - hit_location[0]) > 0.1 ||  fabs(ysum/nhits - hit_location[1]) > 0.1 ||  fabs(zsum/nhits - hit_location[2]) > 0.1)
	      {
		cout << "ALERT! in layer " << layer   << " cluster (x,y,z) and hit (x,y,z) are different!" << endl;
		cout << "     From  geometry object: hit x " << hit_location[0] << " hit y " << hit_location[1] << " hit z " << hit_location[2]  << endl;
		cout << "     From cluster:  nhits " << nhits << " clusx  = " << xsum/nhits << " clusy " << ysum/nhits << " clusz "  <<  zsum/nhits << endl;		
	      }
	  }
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
      geom->find_segment_center(ladder_z_index,
				ladder_phi_index,
				ladder_location);
      double ladderphi = atan2( ladder_location[1], ladder_location[0] );
      ladderphi += geom->get_strip_phi_tilt();
      
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

      const float corr_factor = 1.0; // ladder

      TMatrixF ERR(3,3);
      ERR[0][0] = pow(thickness*invsqrt12*corr_factor,2);
      ERR[0][1] = 0.0;
      ERR[0][2] = 0.0;
      ERR[1][0] = 0.0;
      ERR[1][1] = pow(phisize*invsqrt12*corr_factor,2);
      ERR[1][2] = 0.0;
      ERR[2][0] = 0.0;
      ERR[2][1] = 0.0;
      ERR[2][2] = pow(zsize*invsqrt12*corr_factor,2);
      
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
      R = ROT;
      
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
      
      if (clus_energy > get_threshold_by_layer(layer)) {
	SvtxCluster* ptr = _clusterlist->insert(&clus);
	if (!ptr->isValid()) {
	  static bool first = true;
	  if (first) {
	    cout << PHWHERE << "ERROR: Invalid SvtxClusters are being produced" << endl;
	    ptr->identify();
	    first = false;
	  }
	}

	if(Verbosity() > 1)
	  {
	    // fairly complete sanity check:
	    // Get the list of g4hit positions for this cluster and compare positions	    
	    cout << " For cluster " << ptr->get_id() << " in layer " << ptr->get_layer() << " using e_weighting " << _make_e_weights[layer] << endl;
	    cout << " cluster position: x " << ptr->get_x() << " y " << ptr->get_y() << " z " << ptr->get_z() << endl;
	    cout << " list of SvtxHit id's: " << endl;
	    for (SvtxCluster::HitIter iter = ptr->begin_hits(); iter != ptr->end_hits(); ++iter) {
	      cout << "  " << *iter << " ";
	      SvtxHit *hit = _hits->get(*iter);
	      cout << " cell id from hit = " << hit->get_cellid() << " : " << endl; 
	      PHG4Cell *cell = cells->findCell(hit->get_cellid());
	      double hit_location[3] = {0.0,0.0,0.0};
	      geom->find_strip_center(cell->get_ladder_z_index(),
				      cell->get_ladder_phi_index(),
				      cell->get_zbin(),
				      cell->get_phibin(),
				      hit_location);
	      cout  << "      cell data: "
		    << " layer " << cell->get_layer()
		    << " ladder z index " << cell->get_ladder_z_index()
		    << " ladder phi index " << cell->get_ladder_phi_index()
		    << " ladder z bin " << cell->get_zbin()
		    << " ladder phi bin " << cell->get_phibin()
		    << " strip x,y,z " << hit_location[0] << "  " << hit_location[1] << "  " << hit_location[2] 
		    << " cell edep " << cell->get_edep()
		    << endl;
	      
	      for (PHG4Cell::EdepConstIterator g4iter = cell->get_g4hits().first;
		   g4iter != cell->get_g4hits().second;
		   ++g4iter) {		
		PHG4Hit *g4hit = g4hits->findHit(g4iter->first);
		cout << "      g4hit entry position: x " << g4hit->get_x(0) << " y " << g4hit->get_y(0) << " z " << g4hit->get_z(0) << " edep " << g4hit->get_edep() << " layer " << g4hit->get_layer() << endl;		
		cout << "      g4hit exit position: x " << g4hit->get_x(1) << " y " << g4hit->get_y(1) << " z " << g4hit->get_z(1) << " edep " << g4hit->get_edep() << " layer " << g4hit->get_layer() << endl;		

		// test that there is not a large difference between the cluster position and the entry position(s) of the g4 hits that contributed to it
		if( fabs(ptr->get_x() - g4hit->get_x(0)) > 0.1 ||  fabs(ptr->get_y() - g4hit->get_y(0)) > 0.1 ||  fabs(ptr->get_z() - g4hit->get_z(0)) > 2.0)
		  cout << "        ClusterLadderCells: ALERT! g4hit entry point and cluster (x,y,z) do not agree " << endl;
	      }
	    }
	    cout << endl;
	  }
	
	if (Verbosity()>1) {
	  double radius = sqrt(clusx*clusx+clusy*clusy);
	  double clusphi = atan2(clusy,clusx);
	  cout << "INTT ladder cluster r=" << radius << " phi=" << clusphi << " z=" << clusz << endl;
	  cout << "pos=(" << clus.get_position(0) << ", " << clus.get_position(1)
	       << ", " << clus.get_position(2) << ")" << endl;
	  cout << endl;
	}
      }	else if (Verbosity()>1) {
	double radius = sqrt(clusx*clusx+clusy*clusy);
	double clusphi = atan2(clusy,clusx);
	cout << "removed r=" << radius << " phi=" << clusphi << " z=" << clusz << endl;
	cout << "pos=(" << clus.get_position(0) << ", " << clus.get_position(1)
	     << ", " << clus.get_position(2) << ")" << endl;
	cout << endl;
      } 
    }
  }
  
  return;
}

void PHG4SvtxClusterizer::ClusterMapsLadderCells(PHCompositeNode *topNode) {

  if(Verbosity() > 0)
    cout << "Entering PHG4SvtxClusterizer::ClusterMapsLadderCells " << endl;

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

    if(Verbosity() >4)
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
    sort(cell_list.begin(), cell_list.end(), PHG4SvtxClusterizer::maps_ladder_lessthan);

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

      if(Verbosity() > 2)
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

	if(Verbosity() > 2)
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

	//if(Verbosity() > 2)	
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
	
	if(Verbosity() > 2)
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
      
      if (clus_energy > get_threshold_by_layer(layer)) {
	SvtxCluster* ptr = _clusterlist->insert(&clus);
	if (!ptr->isValid()) {
	  static bool first = true;
	  if (first) {
	    cout << PHWHERE << "ERROR: Invalid SvtxClusters are being produced" << endl;
	    ptr->identify();
	    first = false;
	  }
	}
	if(Verbosity() > 1)
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

	if (Verbosity()>1) {
	  double radius = sqrt(clusx*clusx+clusy*clusy);
	  double clusphi = atan2(clusy,clusx);
	  cout << "clus energy " << clus_energy << " clus_adc " << clus_adc << " r=" << radius << " phi=" << clusphi << " z=" << clusz << endl;
	  cout << "pos=(" << clus.get_position(0) << ", " << clus.get_position(1)
	       << ", " << clus.get_position(2) << ")" << endl;
	  cout << endl;
	}
      }	else if (Verbosity()>1) {
	double radius = sqrt(clusx*clusx+clusy*clusy);
	double clusphi = atan2(clusy,clusx);
	cout << "MAPS ladder cell: removed, clus_energy = " << clus_energy << " below threshold of " <<  get_threshold_by_layer(layer)  << " clus_adc " << clus_adc <<  " r=" << radius << " phi=" << clusphi << " z=" << clusz << endl;
	cout << "pos=(" << clus.get_position(0) << ", " << clus.get_position(1)
	     << ", " << clus.get_position(2) << ")" << endl;
	cout << endl;
      } 
    }
  }
  
  return;
}


void PHG4SvtxClusterizer::PrintClusters(PHCompositeNode *topNode) {

  if (Verbosity() >= 1) {

    SvtxClusterMap *clusterlist = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
    if (!clusterlist) return;
    
    cout << "================= PHG4SvtxClusterizer::process_event() ====================" << endl;
  

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
