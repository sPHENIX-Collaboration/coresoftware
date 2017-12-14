#include "MvtxClusterizer.h"

#include <tracker/TrackerDefs.h>
#include <tracker/TrackerCluster.h>
#include <tracker/TrackerClusterv1.h>
#include <tracker/TrackerClusterContainer.h>
#include <tracker/TrackerHit.h>
#include <tracker/TrackerHitContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeom_MAPS.h>
#include <g4detectors/PHG4CylinderGeom_Siladders.h>
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

static const float twopi = 2.0 * M_PI;

bool MvtxClusterizer::maps_ladder_lessthan(const TrackerHit* lhs,
    const TrackerHit* rhs)
{
  unsigned int lrow = TrackerDefs::MVTXBinning::get_row(lhs->get_hitid());
  unsigned int rrow = TrackerDefs::MVTXBinning::get_row(rhs->get_hitid());

  if ( lrow < rrow ) return true;
  else if ( lrow == rrow )
  {
    unsigned int lcol = TrackerDefs::MVTXBinning::get_col(lhs->get_hitid());
    unsigned int rcol = TrackerDefs::MVTXBinning::get_col(rhs->get_hitid());

    if ( lcol < rcol ) return true;
  }


  return false;
}

bool MvtxClusterizer::maps_ladder_are_adjacent(const TrackerHit* lhs,
    const TrackerHit* rhs)
{
  TrackerDefs::keytype lkey = lhs->get_hitid();
  TrackerDefs::keytype rkey = rhs->get_hitid();

  // want to cluster only within a chip
  if (TrackerDefs::get_layer(lkey) != TrackerDefs::get_layer(rkey))
    return false;
  if (TrackerDefs::MVTXBinning::get_ladder(lkey) != TrackerDefs::MVTXBinning::get_ladder(rkey))
    return false;
  if (TrackerDefs::MVTXBinning::get_chip(lkey) != TrackerDefs::MVTXBinning::get_chip(rkey))
    return false;

  unsigned int lrow = TrackerDefs::MVTXBinning::get_row(lhs->get_hitid());
  unsigned int rrow = TrackerDefs::MVTXBinning::get_row(rhs->get_hitid());
  unsigned int lcol = TrackerDefs::MVTXBinning::get_col(lhs->get_hitid());
  unsigned int rcol = TrackerDefs::MVTXBinning::get_col(rhs->get_hitid());


  if (get_z_clustering(TrackerDefs::get_layer(lkey)))
  {
    if ( fabs(lcol - rcol) <= 1 && fabs(lrow - rrow) <= 1)
      return true;
  }
  else
  {
    if ( fabs(lcol - rcol) == 0 && fabs(lrow - rrow) <= 1)
      return true;
  }

  return false;
}

MvtxClusterizer::MvtxClusterizer(const string &name,
                                 unsigned int min_layer,
                                 unsigned int max_layer) :
  SubsysReco(name),
  _hits(NULL),
  _clusters(NULL),
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
  _hits = findNode::getClass<TrackerHitContainer>(topNode, "TrackerHitContainer");
  if (!_hits) {
    cout << PHWHERE << "ERROR: Can't find node TrackerHitContainer" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  //-----------------
  // Add Cluster Node
  //-----------------

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode
    = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter_dst(dstNode);

  // Create the SVX node if required
  PHCompositeNode* svxNode
    = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode", "SVTX"));
  if (!svxNode) {
    svxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svxNode);
  }

  // Create the Cluster node if required
  TrackerClusterContainer *trackerclusters
    = findNode::getClass<TrackerClusterContainer>(dstNode, "TrackerClusterContainer");
  if (!trackerclusters)
  {
    trackerclusters = new TrackerClusterContainer();
    PHIODataNode<PHObject> *TrackerClusterContainerNode =
      new PHIODataNode<PHObject>(trackerclusters, "TrackerClusterContainer", "PHObject");
    svxNode->addNode(TrackerClusterContainerNode);
  }

  //----------------
  // Report Settings
  //----------------

  if (verbosity > 0) {
    cout << "====================== MvtxClusterizer::InitRun() =====================" << endl;
    cout << " Fraction of expected thickness MIP energy = " << _fraction_of_mip << endl;
    for (std::map<int, float>::iterator iter = _thresholds_by_layer.begin();
         iter != _thresholds_by_layer.end();
         ++iter) {
      cout << " Cluster Threshold in Layer #" << iter->first << " = " << 1.0e6 * iter->second << " keV" << endl;
    }
    for (std::map<int, bool>::iterator iter = _make_z_clustering.begin();
         iter != _make_z_clustering.end();
         ++iter) {
      cout << " Z-dimension Clustering in Layer #" << iter->first << " = " << boolalpha << iter->second << noboolalpha << endl;
    }
    for (std::map<int, bool>::iterator iter = _make_e_weights.begin();
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

  _clusters = findNode::getClass<TrackerClusterContainer>(topNode, "TrackerClusterContainer");
  if (!_clusters)
  {
    cout << PHWHERE << " ERROR: Can't find TrackerClusterContainer." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  _clusters->Reset();

  _hits = findNode::getClass<TrackerHitContainer>(topNode, "TrackerHitContainer");
  if (!_hits)
  {
    cout << PHWHERE << " ERROR: Can't find TrackerHitContainer." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  ClusterMapsLadderHits(topNode);

  PrintClusters(topNode);

  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

void MvtxClusterizer::ClusterMapsLadderHits(PHCompositeNode *topNode)
{

  if (verbosity > 0)
    cout << "Entering MvtxClusterizer::ClusterMapsLadderCells " << endl;

  //----------
  // Get Nodes
  //----------

  // get the SVX geometry object
  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MAPS");
  if (!geom_container) return;

  //-----------
  // Clustering
  //-----------


  // Loop over each layer, ladder, chip and cluster hits only within the chip
  // D. McGlinchey -- Loop over only layer for now, can make this more
  //                  efficient once we change geometry
  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for (PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {

    int layer = layeriter->second->get_layer();

    if ((unsigned int)layer < _min_layer) continue;
    if ((unsigned int)layer > _max_layer) continue;

    // get the layer geometry helper
    PHG4CylinderGeom_MAPS *geom = (PHG4CylinderGeom_MAPS*) geom_container->GetLayerGeom(layer);

    // get all the hits in this layer and put them into a vector
    vector<TrackerHit*> hit_list;
    TrackerHitContainer::ConstRange mvtxhitrange = _hits->getHits(TrackerDefs::TRACKERID::mvtx_id, layer);
    for (TrackerHitContainer::ConstIterator iter = mvtxhitrange.first;
         iter != mvtxhitrange.second;
         ++iter)
    {
      // D. McGlinchey - Figure out the correct const'ness to use here ...
      // hit_list.insert(iter->second);
      hit_list.push_back(_hits->findHit(iter->second->get_hitid()));
    }

    // i'm not sure this sorting is ever really used
    sort(hit_list.begin(), hit_list.end(), MvtxClusterizer::maps_ladder_lessthan);

    typedef adjacency_list <vecS, vecS, undirectedS> Graph;
    Graph G;

    for (unsigned int i = 0; i < hit_list.size(); i++)
    {
      for (unsigned int j = i + 1; j < hit_list.size(); j++)
      {
        if (maps_ladder_are_adjacent(hit_list[i], hit_list[j]) )
          add_edge(i, j, G);
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
    set<int> cluster_ids; // unique components
    multimap<int, TrackerHit*> clusters;
    for (unsigned int i = 0; i < component.size(); i++)
    {
      cluster_ids.insert( component[i] );
      clusters.insert( make_pair(component[i], hit_list[i]) );
    }

    //
    for (set<int>::iterator clusiter = cluster_ids.begin();
         clusiter != cluster_ids.end();
         clusiter++ )
    {

      int clusid = *clusiter;
      pair<multimap<int, TrackerHit*>::iterator,
           multimap<int, TrackerHit*>::iterator> clusrange = clusters.equal_range(clusid);

      multimap<int, TrackerHit*>::iterator mapiter = clusrange.first;

      if (verbosity > 2)
        cout << "Filling cluster id " << clusid << " in  layer " << layer << endl;

      TrackerClusterv1* clus;

      // determine the size of the cluster in phi and z
      // useful for track fitting the cluster
      set<unsigned int> phibins;
      set<unsigned int> zbins;
      for (mapiter = clusrange.first; mapiter != clusrange.second; mapiter++ )
      {
        TrackerHit* hit = mapiter->second;

        unsigned int binphi = TrackerDefs::MVTXBinning::get_row(hit->get_hitid());
        phibins.insert(binphi);
        unsigned int binz = TrackerDefs::MVTXBinning::get_col(hit->get_hitid());
        zbins.insert(binz);
      }

      float thickness = geom->get_pixel_thickness();
      float pitch = geom->get_pixel_x();
      float length = geom->get_pixel_z();
      float phisize = phibins.size() * pitch;
      float zsize = zbins.size() * length;
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

      for (mapiter = clusrange.first; mapiter != clusrange.second; mapiter++ )
      {
        TrackerHit* hit = mapiter->second;

        clus->insert_hit(hit->get_hitid());

        // find the center of the pixel in world coords
        TrackerDefs::keytype hitid = hit->get_hitid();
        unsigned int row = TrackerDefs::MVTXBinning::get_row(hitid);
        unsigned int col = TrackerDefs::MVTXBinning::get_col(hitid);
        int pixel_number = geom->get_pixel_number_from_xbin_zbin(row, col);

        // These will be used later to get the sensor position so that the sensor phi can be calculated
        stave_index = TrackerDefs::MVTXBinning::get_ladder(hitid);
        half_stave_index = 0; //cell->get_half_stave_index();
        module_index = 0; //cell->get_module_index();
        chip_index = TrackerDefs::MVTXBinning::get_chip(hitid);

        TVector3 local_coords = geom->get_local_coords_from_pixel(pixel_number);
        TVector3 world_coords =
          geom->get_world_from_local_coords(stave_index, half_stave_index, module_index, chip_index, local_coords);
        double hit_location[3] = {world_coords.X(), world_coords.Y(), world_coords.Z()};


        xsum += hit_location[0];
        ysum += hit_location[1];
        zsum += hit_location[2];

        ++nhits;
      }

      double clusx = xsum / nhits;
      double clusy = ysum / nhits;
      double clusz = zsum / nhits;

      double ladder_location[3] = {0.0, 0.0, 0.0};
      // returns the center of the sensor in world coordinates - used to get the ladder phi location
      geom->find_sensor_center(stave_index, half_stave_index, module_index, chip_index, ladder_location);
      double ladderphi = atan2( ladder_location[1], ladder_location[0] );
      ladderphi += geom->get_stave_phi_tilt();

      //cout << "sensor center = " << ladder_location[0] << " " << ladder_location[1] << " " << ladder_location[2] << endl;

      clus->set_position(0, clusx);
      clus->set_position(1, clusy);
      clus->set_position(2, clusz);

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

      clus->set_size( 0 , 0 , COVAR_DIM[0][0] );
      clus->set_size( 0 , 1 , COVAR_DIM[0][1] );
      clus->set_size( 0 , 2 , COVAR_DIM[0][2] );
      clus->set_size( 1 , 0 , COVAR_DIM[1][0] );
      clus->set_size( 1 , 1 , COVAR_DIM[1][1] );
      clus->set_size( 1 , 2 , COVAR_DIM[1][2] );
      clus->set_size( 2 , 0 , COVAR_DIM[2][0] );
      clus->set_size( 2 , 1 , COVAR_DIM[2][1] );
      clus->set_size( 2 , 2 , COVAR_DIM[2][2] );

      TMatrixF COVAR_ERR(3, 3);
      COVAR_ERR = R * ERR * R_T;

      clus->set_error( 0 , 0 , COVAR_ERR[0][0] );
      clus->set_error( 0 , 1 , COVAR_ERR[0][1] );
      clus->set_error( 0 , 2 , COVAR_ERR[0][2] );
      clus->set_error( 1 , 0 , COVAR_ERR[1][0] );
      clus->set_error( 1 , 1 , COVAR_ERR[1][1] );
      clus->set_error( 1 , 2 , COVAR_ERR[1][2] );
      clus->set_error( 2 , 0 , COVAR_ERR[2][0] );
      clus->set_error( 2 , 1 , COVAR_ERR[2][1] );
      clus->set_error( 2 , 2 , COVAR_ERR[2][2] );


      // make the unique cluster key
      TrackerDefs::keytype cid =
        TrackerDefs::MVTXBinning::gencluskey(
          TrackerDefs::TRACKERID::mvtx_id,
          layer, stave_index, chip_index,
          clusid);
      clus->set_id(cid);

      _clusters->AddClusterSpecifyKey(cid, clus);

    }
  } // layeriter





  return;
}


void MvtxClusterizer::PrintClusters(PHCompositeNode * topNode) {

  if (verbosity >= 1) 
  {

    TrackerClusterContainer::ConstRange clusrange = _clusters->getClusters(TrackerDefs::TRACKERID::mvtx_id);

    cout << "================= MvtxClusterizer::process_event() ====================" << endl;


    cout << " Found and recorded the following clusters: " << endl;

    unsigned int iclus = 0;
    for (TrackerClusterContainer::ConstIterator iter = clusrange.first;
         iter != clusrange.second;
         ++iter)
    {

      TrackerCluster* clus = iter->second;
      cout << iclus << endl;
      clus->identify();
      ++iclus;
    }

    cout << "===========================================================================" << endl;
  }

  return;
}
