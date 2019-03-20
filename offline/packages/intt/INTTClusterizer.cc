#include "INTTClusterizer.h"
#include "CylinderGeomINTT.h"

#include <trackbase_historic/SvtxCluster.h>
#include <trackbase_historic/SvtxClusterMap.h>
#include <trackbase_historic/SvtxClusterMap_v1.h>
#include <trackbase_historic/SvtxCluster_v1.h>
#include <trackbase_historic/SvtxHit.h>
#include <trackbase_historic/SvtxHitMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>

#include <TMatrixF.h>
#include <TVector3.h>

#define BOOST_NO_HASH  // Our version of boost.graph is incompatible with GCC-4.3 w/o this flag
#include <boost/bind.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
using namespace boost;

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace std;

static const float twopi = 2.0 * M_PI;

bool INTTClusterizer::lessthan(const PHG4Cell* lhs,
                               const PHG4Cell* rhs)
{
  int lhsphibin = PHG4CellDefs::SizeBinning::get_phibin(lhs->get_cellid());
  int rhsphibin = PHG4CellDefs::SizeBinning::get_phibin(rhs->get_cellid());
  int lhszbin = PHG4CellDefs::SizeBinning::get_zbin(lhs->get_cellid());
  int rhszbin = PHG4CellDefs::SizeBinning::get_zbin(rhs->get_cellid());

  if (lhsphibin < rhsphibin)
    return true;
  else if (lhsphibin == rhsphibin)
  {
    if (lhszbin < rhszbin) return true;
  }

  return false;
}

bool INTTClusterizer::ladder_lessthan(const PHG4Cell* lhs,
                                      const PHG4Cell* rhs)
{
  if (lhs->get_ladder_z_index() == rhs->get_ladder_z_index() &&
      lhs->get_ladder_phi_index() == rhs->get_ladder_phi_index())
  {
    if (lhs->get_phibin() < rhs->get_phibin())
      return true;
    else if (lhs->get_phibin() == rhs->get_phibin())
    {
      if (lhs->get_zbin() < rhs->get_zbin()) return true;
    }
  }
  else
  {
    if (lhs->get_zbin() < rhs->get_zbin()) return true;
  }

  return false;
}

bool INTTClusterizer::are_adjacent(const PHG4Cell* lhs,
                                   const PHG4Cell* rhs,
                                   const int& nphibins)
{
  int lhs_layer = lhs->get_layer();
  int rhs_layer = rhs->get_layer();
  if (lhs_layer != rhs_layer) return false;

  int lhsphibin = PHG4CellDefs::SizeBinning::get_phibin(lhs->get_cellid());
  int rhsphibin = PHG4CellDefs::SizeBinning::get_phibin(rhs->get_cellid());
  int lhszbin = PHG4CellDefs::SizeBinning::get_zbin(lhs->get_cellid());
  int rhszbin = PHG4CellDefs::SizeBinning::get_zbin(rhs->get_cellid());
  if (get_z_clustering(lhs_layer))
  {
    if (fabs(lhszbin - rhszbin) <= 1)
    {
      if (fabs(lhsphibin - rhsphibin) <= 1)
      {
        return true;
      }
      else if (lhsphibin == 0 || rhsphibin == 0)
      {
        if (fabs(lhsphibin - rhsphibin) == (nphibins - 1))
          return true;
      }
    }
  }
  else
  {
    if (fabs(lhszbin - rhszbin) == 0)
    {
      if (fabs(lhsphibin - rhsphibin) <= 1)
      {
        return true;
      }
      else if (lhsphibin == 0 || rhsphibin == 0)
      {
        if (fabs(lhsphibin - rhsphibin) == (nphibins - 1))
          return true;
      }
    }
  }

  return false;
}

bool INTTClusterizer::ladder_are_adjacent(const PHG4Cell* lhs, const PHG4Cell* rhs)
{
  if (Verbosity() > 2)
  {
    cout << " lhs layer " << lhs->get_layer()
         << " lhs ladder z index " << lhs->get_ladder_z_index()
         << " lhs ladder phi index " << lhs->get_ladder_phi_index()
         << " lhs z bin " << lhs->get_zbin()
         << " lhs phi bin " << lhs->get_phibin()
         << endl;

    cout << " rhs layer " << rhs->get_layer()
         << " rhs ladder z index " << rhs->get_ladder_z_index()
         << " rhs ladder phi index " << rhs->get_ladder_phi_index()
         << " rhs z bin " << rhs->get_zbin()
         << " rhs phi bin " << rhs->get_phibin()
         << endl;
  }

  int lhs_layer = lhs->get_layer();
  int rhs_layer = rhs->get_layer();
  if (lhs_layer != rhs_layer) return false;

  if (!(lhs->get_ladder_z_index() == rhs->get_ladder_z_index() &&
        lhs->get_ladder_phi_index() == rhs->get_ladder_phi_index())) return false;

  if (get_z_clustering(lhs_layer))
  {
    if (fabs(lhs->get_zbin() - rhs->get_zbin()) <= 1)
    {
      if (fabs(lhs->get_phibin() - rhs->get_phibin()) <= 1)
      {
        return true;
      }
    }
  }
  else
  {
    if (fabs(lhs->get_zbin() - rhs->get_zbin()) == 0)
    {
      if (fabs(lhs->get_phibin() - rhs->get_phibin()) <= 1)
      {
        if (Verbosity() > 2) cout << "    accepted " << endl;
        return true;
      }
    }
  }

  return false;
}

INTTClusterizer::INTTClusterizer(const string& name,
                                 unsigned int min_layer,
                                 unsigned int max_layer)
  : SubsysReco(name)
  , _hits(NULL)
  , _clusterlist(NULL)
  , _fraction_of_mip(0.5)
  , _thresholds_by_layer()
  , _make_z_clustering()
  , _make_e_weights()
  , _min_layer(min_layer)
  , _max_layer(max_layer)
  , _timer(PHTimeServer::get()->insert_new(name))
{
}

int INTTClusterizer::InitRun(PHCompositeNode* topNode)
{
  // get node containing the digitized hits
  _hits = findNode::getClass<SvtxHitMap>(topNode, "SvtxHitMap");
  if (!_hits)
  {
    cout << PHWHERE << "ERROR: Can't find node SvtxHitMap" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

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

  // Create the SVX node if required
  PHCompositeNode* svxNode = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode", "SVTX"));
  if (!svxNode)
  {
    svxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svxNode);
  }

  // Create the Cluster node if required
  SvtxClusterMap* svxclusters = findNode::getClass<SvtxClusterMap>(dstNode, "SvtxClusterMap");
  if (!svxclusters)
  {
    svxclusters = new SvtxClusterMap_v1();
    PHIODataNode<PHObject>* SvtxClusterMapNode =
        new PHIODataNode<PHObject>(svxclusters, "SvtxClusterMap", "PHObject");
    svxNode->addNode(SvtxClusterMapNode);
  }

  //---------------------
  // Calculate Thresholds
  //---------------------

  CalculateLadderThresholds(topNode);

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    cout << "====================== INTTClusterizer::InitRun() =====================" << endl;
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

int INTTClusterizer::process_event(PHCompositeNode* topNode)
{
  _timer.get()->restart();

  _clusterlist = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
  if (!_clusterlist)
  {
    cout << PHWHERE << " ERROR: Can't find SvtxClusterMap." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  ClusterLadderCells(topNode);

  PrintClusters(topNode);

  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

void INTTClusterizer::CalculateLadderThresholds(PHCompositeNode* topNode)
{
  PHG4CellContainer* cells = findNode::getClass<PHG4CellContainer>(topNode, "G4CELL_INTT");
  if (!cells) return;

  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if (!geom_container) return;

  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for (PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    // index mapping
    int layer = layeriter->second->get_layer();

    // cluster threshold
    float thickness = (layeriter->second)->get_thickness();
    float threshold = _fraction_of_mip * 0.003876 * thickness;
    _thresholds_by_layer.insert(std::make_pair(layer, threshold));

    // fill in a default z_clustering value if not present
    if (_make_z_clustering.find(layer) == _make_z_clustering.end())
    {
      _make_z_clustering.insert(std::make_pair(layer, true));
    }

    if (_make_e_weights.find(layer) == _make_e_weights.end())
    {
      _make_e_weights.insert(std::make_pair(layer, false));
    }
  }

  return;
}

void INTTClusterizer::ClusterLadderCells(PHCompositeNode* topNode)
{
  if (Verbosity() > 0)
    cout << "Entering INTTClusterizer::ClusterLadderCells " << endl;

  //----------
  // Get Nodes
  //----------

  // get the SVX geometry object
  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if (!geom_container) return;

  PHG4HitContainer* g4hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  if (!g4hits) return;

  PHG4CellContainer* cells = findNode::getClass<PHG4CellContainer>(topNode, "G4CELL_INTT");
  if (!cells) return;

  //-----------
  // Clustering
  //-----------

  // sort hits layer by layer
  std::multimap<int, SvtxHit*> layer_hits_mmap;
  for (SvtxHitMap::Iter iter = _hits->begin();
       iter != _hits->end();
       ++iter)
  {
    SvtxHit* hit = iter->second;
    layer_hits_mmap.insert(make_pair(hit->get_layer(), hit));
  }

  // this does nothing!
  for (std::multimap<int, SvtxHit*>::iterator it = layer_hits_mmap.begin(); it != layer_hits_mmap.end(); ++it)
  {
    if ((*it).first < (int) _min_layer || (*it).first > (int) _max_layer)
      continue;
  }

  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for (PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    int layer = layeriter->second->get_layer();
    if (Verbosity() > 0)
      cout << " layer loop, current layer = " << layer << endl;

    if ((unsigned int) layer < _min_layer) continue;
    if ((unsigned int) layer > _max_layer) continue;

    std::map<PHG4Cell*, SvtxHit*> cell_hit_map;
    vector<PHG4Cell*> cell_list;
    for (std::multimap<int, SvtxHit*>::iterator hiter = layer_hits_mmap.lower_bound(layer);
         hiter != layer_hits_mmap.upper_bound(layer);
         ++hiter)
    {
      SvtxHit* hit = hiter->second;
      PHG4Cell* cell = cells->findCell(hit->get_cellid());
      if (Verbosity() > 2)
      {
        cout << "adding cell to cell_hit_map: ";
        cell->print();
      }
      cell_list.push_back(cell);
      cell_hit_map.insert(make_pair(cell, hit));
    }

    if (cell_list.size() == 0) continue;  // if no cells, go to the next layer

    // i'm not sure this sorting is ever really used
    sort(cell_list.begin(), cell_list.end(), INTTClusterizer::ladder_lessthan);

    typedef adjacency_list<vecS, vecS, undirectedS> Graph;
    Graph G;

    // Find adjacent cell
    if (Verbosity() > 2) cout << "Find adjacent cells for layer " << layer << endl;
    for (unsigned int i = 0; i < cell_list.size(); i++)
    {
      for (unsigned int j = i + 1; j < cell_list.size(); j++)
      {
        if (Verbosity() > 2)
        {
          cout << "compare cells " << i << " and " << j << endl;
          cell_list[i]->print();
          cell_list[j]->print();
        }
        if (ladder_are_adjacent(cell_list[i], cell_list[j]))
        {
          add_edge(i, j, G);
          if (Verbosity() > 2)
          {
            cout << "Found edge " << i << "   " << j << endl;
            cell_list[i]->print();
            cell_list[j]->print();
          }
        }
      }

      add_edge(i, i, G);
    }
    if (Verbosity() > 2) cout << "finished looking for adjacent cells for layer " << layer << endl;

    // Find the connections between the vertices of the graph (vertices are the rawhits,
    // connections are made when they are adjacent to one another)
    vector<int> component(num_vertices(G));

    // this is the actual clustering, performed by boost
    connected_components(G, &component[0]);

    // Loop over the components(hit cells) compiling a list of the
    // unique connected groups (ie. clusters).
    set<int> cluster_ids;  // unique components
    multimap<int, PHG4Cell*> clusters;
    for (unsigned int i = 0; i < component.size(); i++)
    {
      cluster_ids.insert(component[i]);
      clusters.insert(make_pair(component[i], cell_list[i]));
    }

    //
    for (set<int>::iterator clusiter = cluster_ids.begin();
         clusiter != cluster_ids.end();
         ++clusiter)
    {
      int clusid = *clusiter;
      pair<multimap<int, PHG4Cell*>::iterator,
           multimap<int, PHG4Cell*>::iterator>
          clusrange = clusters.equal_range(clusid);

      multimap<int, PHG4Cell*>::iterator mapiter = clusrange.first;

      int layer = mapiter->second->get_layer();
      //PHG4CylinderGeom* geom = geom_container->GetLayerGeom(layer);
      CylinderGeomINTT* geom = dynamic_cast<CylinderGeomINTT*>(geom_container->GetLayerGeom(layer));

      SvtxCluster_v1 clus;
      clus.set_layer(layer);
      float clus_energy = 0.0;
      unsigned int clus_adc = 0;

      // determine the size of the cluster in phi and z
      // useful for track fitting the cluster

      set<int> phibins;
      set<int> zbins;
      for (mapiter = clusrange.first; mapiter != clusrange.second; ++mapiter)
      {
        PHG4Cell* cell = mapiter->second;

        phibins.insert(cell->get_phibin());
        zbins.insert(cell->get_zbin());
      }

      float thickness = geom->get_thickness();
      float pitch = geom->get_strip_y_spacing();
      float length = geom->get_strip_z_spacing();
      float phisize = phibins.size() * pitch;
      float zsize = zbins.size() * length;

      // tilt refers to a rotation around the radial vector from the origin, and this is zero for the INTT ladders
      float tilt = 0;  //geom->get_strip_tilt();

      // determine the cluster position...
      double xsum = 0.0;
      double ysum = 0.0;
      double zsum = 0.0;
      unsigned nhits = 0;

      int ladder_z_index = -1;
      int ladder_phi_index = -1;

      for (mapiter = clusrange.first; mapiter != clusrange.second; ++mapiter)
      {
        PHG4Cell* cell = mapiter->second;
        SvtxHit* hit = cell_hit_map[cell];
        if (Verbosity() > 0)
        {
          cout << "Add hit: ";
          hit->identify();
          cout << " cell is ";
          cell->print();
        }
        clus.insert_hit(hit->get_id());

        clus_energy += hit->get_e();
        clus_adc += hit->get_adc();

        double hit_location[3] = {0.0, 0.0, 0.0};
        geom->find_strip_center(cell->get_ladder_z_index(),
                                cell->get_ladder_phi_index(),
                                cell->get_zbin(),
                                cell->get_phibin(),
                                hit_location);

        ladder_z_index = cell->get_ladder_z_index();
        ladder_phi_index = cell->get_ladder_phi_index();

        if (_make_e_weights[layer])
        {
          xsum += hit_location[0] * hit->get_e();
          ysum += hit_location[1] * hit->get_e();
          zsum += hit_location[2] * hit->get_e();
        }
        else
        {
          xsum += hit_location[0];
          ysum += hit_location[1];
          zsum += hit_location[2];
        }
        ++nhits;
        if (Verbosity() > 2) cout << "     nhits = " << nhits << endl;
        if (Verbosity() > 2)
        {
          cout << "  From  geometry object: hit x " << hit_location[0] << " hit y " << hit_location[1] << " hit z " << hit_location[2] << endl;
          cout << "     nhits " << nhits << " clusx  = " << xsum / nhits << " clusy " << ysum / nhits << " clusz " << zsum / nhits << endl;

          if (fabs(xsum / nhits - hit_location[0]) > 0.1 || fabs(ysum / nhits - hit_location[1]) > 0.1 || fabs(zsum / nhits - hit_location[2]) > 0.1)
          {
            cout << "ALERT! in layer " << layer << " cluster (x,y,z) and hit (x,y,z) are different!" << endl;
            cout << "     From  geometry object: hit x " << hit_location[0] << " hit y " << hit_location[1] << " hit z " << hit_location[2] << endl;
            cout << "     From cluster:  nhits " << nhits << " clusx  = " << xsum / nhits << " clusy " << ysum / nhits << " clusz " << zsum / nhits << endl;
          }
        }
      }

      double clusx = NAN;
      double clusy = NAN;
      double clusz = NAN;

      if (_make_e_weights[layer])
      {
        clusx = xsum / clus_energy;
        clusy = ysum / clus_energy;
        clusz = zsum / clus_energy;
      }
      else
      {
        clusx = xsum / nhits;
        clusy = ysum / nhits;
        clusz = zsum / nhits;
      }

      double ladder_location[3] = {0.0, 0.0, 0.0};
      geom->find_segment_center(ladder_z_index,
                                ladder_phi_index,
                                ladder_location);
      double ladderphi = atan2(ladder_location[1], ladder_location[0]);
      ladderphi += geom->get_strip_phi_tilt();

      clus.set_position(0, clusx);
      clus.set_position(1, clusy);
      clus.set_position(2, clusz);

      clus.set_e(clus_energy);
      clus.set_adc(clus_adc);

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
      R = ROT;

      TMatrixF R_T(3, 3);
      R_T.Transpose(R);

      TMatrixF COVAR_DIM(3, 3);
      COVAR_DIM = R * DIM * R_T;

      clus.set_size(0, 0, COVAR_DIM[0][0]);
      clus.set_size(0, 1, COVAR_DIM[0][1]);
      clus.set_size(0, 2, COVAR_DIM[0][2]);
      clus.set_size(1, 0, COVAR_DIM[1][0]);
      clus.set_size(1, 1, COVAR_DIM[1][1]);
      clus.set_size(1, 2, COVAR_DIM[1][2]);
      clus.set_size(2, 0, COVAR_DIM[2][0]);
      clus.set_size(2, 1, COVAR_DIM[2][1]);
      clus.set_size(2, 2, COVAR_DIM[2][2]);

      TMatrixF COVAR_ERR(3, 3);
      COVAR_ERR = R * ERR * R_T;

      clus.set_error(0, 0, COVAR_ERR[0][0]);
      clus.set_error(0, 1, COVAR_ERR[0][1]);
      clus.set_error(0, 2, COVAR_ERR[0][2]);
      clus.set_error(1, 0, COVAR_ERR[1][0]);
      clus.set_error(1, 1, COVAR_ERR[1][1]);
      clus.set_error(1, 2, COVAR_ERR[1][2]);
      clus.set_error(2, 0, COVAR_ERR[2][0]);
      clus.set_error(2, 1, COVAR_ERR[2][1]);
      clus.set_error(2, 2, COVAR_ERR[2][2]);

      if (clus_energy > get_threshold_by_layer(layer))
      {
        SvtxCluster* ptr = _clusterlist->insert(&clus);
        if (!ptr->isValid())
        {
          static bool first = true;
          if (first)
          {
            cout << PHWHERE << "ERROR: Invalid SvtxClusters are being produced" << endl;
            ptr->identify();
            first = false;
          }
        }

        if (Verbosity() > 1)
        {
          // fairly complete sanity check:
          // Get the list of g4hit positions for this cluster and compare positions
          cout << " For cluster " << ptr->get_id() << " in layer " << ptr->get_layer() << " using e_weighting " << _make_e_weights[layer] << endl;
          cout << " cluster position: x " << ptr->get_x() << " y " << ptr->get_y() << " z " << ptr->get_z() << endl;
          cout << " list of SvtxHit id's: " << endl;
          for (SvtxCluster::HitIter iter = ptr->begin_hits(); iter != ptr->end_hits(); ++iter)
          {
            cout << "  " << *iter << " ";
            SvtxHit* hit = _hits->get(*iter);
            cout << " cell id from hit = " << hit->get_cellid() << " : " << endl;
            PHG4Cell* cell = cells->findCell(hit->get_cellid());
            double hit_location[3] = {0.0, 0.0, 0.0};
            geom->find_strip_center(cell->get_ladder_z_index(),
                                    cell->get_ladder_phi_index(),
                                    cell->get_zbin(),
                                    cell->get_phibin(),
                                    hit_location);
            cout << "      cell data: "
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
                 ++g4iter)
            {
              PHG4Hit* g4hit = g4hits->findHit(g4iter->first);
              cout << "      g4hit entry position: x " << g4hit->get_x(0) << " y " << g4hit->get_y(0) << " z " << g4hit->get_z(0) << " edep " << g4hit->get_edep() << " layer " << g4hit->get_layer() << endl;
              cout << "      g4hit exit position: x " << g4hit->get_x(1) << " y " << g4hit->get_y(1) << " z " << g4hit->get_z(1) << " edep " << g4hit->get_edep() << " layer " << g4hit->get_layer() << endl;

              // test that there is not a large difference between the cluster position and the entry position(s) of the g4 hits that contributed to it
              if (fabs(ptr->get_x() - g4hit->get_x(0)) > 0.1 || fabs(ptr->get_y() - g4hit->get_y(0)) > 0.1 || fabs(ptr->get_z() - g4hit->get_z(0)) > 2.0)
                cout << "        ClusterLadderCells: ALERT! g4hit entry point and cluster (x,y,z) do not agree " << endl;
            }
          }
          cout << endl;
        }

        if (Verbosity() > 1)
        {
          double radius = sqrt(clusx * clusx + clusy * clusy);
          double clusphi = atan2(clusy, clusx);
          cout << "INTT ladder cluster r=" << radius << " phi=" << clusphi << " z=" << clusz << endl;
          cout << "pos=(" << clus.get_position(0) << ", " << clus.get_position(1)
               << ", " << clus.get_position(2) << ")" << endl;
          cout << endl;
        }
      }
      else if (Verbosity() > 1)
      {
        double radius = sqrt(clusx * clusx + clusy * clusy);
        double clusphi = atan2(clusy, clusx);
        cout << "removed r=" << radius << " phi=" << clusphi << " z=" << clusz << endl;
        cout << "pos=(" << clus.get_position(0) << ", " << clus.get_position(1)
             << ", " << clus.get_position(2) << ")" << endl;
        cout << endl;
      }
    }
  }

  return;
}

void INTTClusterizer::PrintClusters(PHCompositeNode* topNode)
{
  if (Verbosity() >= 1)
  {
    SvtxClusterMap* clusterlist = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
    if (!clusterlist) return;

    cout << "================= INTTClusterizer::process_event() ====================" << endl;

    cout << " Found and recorded the following " << clusterlist->size() << " clusters: " << endl;

    unsigned int icluster = 0;
    for (SvtxClusterMap::Iter iter = clusterlist->begin();
         iter != clusterlist->end();
         ++iter)
    {
      SvtxCluster* cluster = iter->second;
      cout << icluster << " of " << clusterlist->size() << endl;
      cluster->identify();
      ++icluster;
    }

    cout << "===========================================================================" << endl;
  }

  return;
}
