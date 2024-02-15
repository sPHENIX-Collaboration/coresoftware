#include "InttClusterizer.h"
#include "CylinderGeomIntt.h"

#include <trackbase/InttDefs.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterCrossingAssocv1.h>
#include <trackbase/TrkrClusterHitAssocv3.h>
#include <trackbase/TrkrClusterv5.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <trackbase/ClusHitsVerbosev1.h>
#include <trackbase/RawHit.h>
#include <trackbase/RawHitSet.h>
#include <trackbase/RawHitSetContainer.h>

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#pragma GCC diagnostic pop

#include <array>
#include <cmath>
#include <iostream>
#include <memory>  // for unique_ptr, make_...
#include <set>
#include <vector>  // for vector

namespace
{

  /// convenience square method
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }
}  // namespace

bool InttClusterizer::ladder_are_adjacent(const std::pair<TrkrDefs::hitkey, TrkrHit*>& lhs, const std::pair<TrkrDefs::hitkey, TrkrHit*>& rhs, const int layer)
{
  if (get_z_clustering(layer))
  {
    if (fabs(InttDefs::getCol(lhs.first) - InttDefs::getCol(rhs.first)) <= 1)
    {
      if (fabs(InttDefs::getRow(lhs.first) - InttDefs::getRow(rhs.first)) <= 1)
      {
        return true;
      }
    }
  }
  else if (fabs(InttDefs::getCol(lhs.first) - InttDefs::getCol(rhs.first)) == 0)
  {
    if (fabs(InttDefs::getRow(lhs.first) - InttDefs::getRow(rhs.first)) <= 1)
    {
      return true;
    }
  }

  return false;
}

bool InttClusterizer::ladder_are_adjacent(RawHit* lhs, RawHit* rhs, const int layer)
{
  if (get_z_clustering(layer))
  {
    if (fabs(lhs->getPhiBin() - rhs->getPhiBin()) <= 1)  // col
    {
      if (fabs(lhs->getTBin() - rhs->getTBin()) <= 1)  // Row
      {
        return true;
      }
    }
  }
  else if (fabs(lhs->getPhiBin() - rhs->getPhiBin()) <= 1)
  {
    if (fabs(lhs->getTBin() - rhs->getTBin()) == 0)
    {
      return true;
    }
  }

  return false;
}

InttClusterizer::InttClusterizer(const std::string& name,
                                 unsigned int /*min_layer*/,
                                 unsigned int /*max_layer*/)
  : SubsysReco(name)
{
}

int InttClusterizer::InitRun(PHCompositeNode* topNode)
{
  /*
  // get node containing the digitized hits
  _hits = findNode::getClass<SvtxHitMap>(topNode, "SvtxHitMap");
  if (!_hits)
  {
    std::cout << PHWHERE << "ERROR: Can't find node SvtxHitMap" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  */

  //-----------------
  // Add Cluster Node
  //-----------------

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter_dst(dstNode);

  // Create the Cluster and association nodes if required
  auto trkrclusters = findNode::getClass<TrkrClusterContainer>(dstNode, "TRKR_CLUSTER");
  if (!trkrclusters)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode* DetNode =
        dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    trkrclusters = new TrkrClusterContainerv4;
    PHIODataNode<PHObject>* TrkrClusterContainerNode =
        new PHIODataNode<PHObject>(trkrclusters, "TRKR_CLUSTER", "PHObject");
    DetNode->addNode(TrkrClusterContainerNode);
  }

  auto clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!clusterhitassoc)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode* DetNode =
        dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    clusterhitassoc = new TrkrClusterHitAssocv3;
    PHIODataNode<PHObject>* newNode = new PHIODataNode<PHObject>(clusterhitassoc, "TRKR_CLUSTERHITASSOC", "PHObject");
    DetNode->addNode(newNode);
  }

  // Add the multimap holding the INTT cluster-crossing associations
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode* DetNode =
        dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    auto clustercrossingassoc = new TrkrClusterCrossingAssocv1;
    PHIODataNode<PHObject>* newNode = new PHIODataNode<PHObject>(clustercrossingassoc, "TRKR_CLUSTERCROSSINGASSOC", "PHObject");
    DetNode->addNode(newNode);
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
    std::cout << "====================== InttClusterizer::InitRun() =====================" << std::endl;
    std::cout << " Fraction of expected thickness MIP energy = " << _fraction_of_mip << std::endl;
    for (auto& iter_tmp : _thresholds_by_layer)
    {
      std::cout << " Cluster Threshold in Layer #" << iter_tmp.first << " = " << 1.0e6 * iter_tmp.second << " keV" << std::endl;
    }
    for (auto& iter_tmp : _make_z_clustering)
    {
      std::cout << " Z-dimension Clustering in Layer #" << iter_tmp.first << " = " << std::boolalpha << iter_tmp.second << std::noboolalpha << std::endl;
    }
    for (auto& _make_e_weight : _make_e_weights)
    {
      std::cout << " Energy weighting clusters in Layer #" << _make_e_weight.first << " = " << std::boolalpha << _make_e_weight.second << std::noboolalpha << std::endl;
    }
    std::cout << "===========================================================================" << std::endl;
  }

  if (record_ClusHitsVerbose)
  {
    // get the node
    mClusHitsVerbose = findNode::getClass<ClusHitsVerbosev1>(topNode, "Trkr_SvtxClusHitsVerbose");
    if (!mClusHitsVerbose)
    {
      PHNodeIterator dstiter(dstNode);
      auto DetNode = dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode", "TRKR"));
      if (!DetNode)
      {
        DetNode = new PHCompositeNode("TRKR");
        dstNode->addNode(DetNode);
      }
      mClusHitsVerbose = new ClusHitsVerbosev1();
      auto newNode = new PHIODataNode<PHObject>(mClusHitsVerbose, "Trkr_SvtxClusHitsVerbose", "PHObject");
      DetNode->addNode(newNode);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttClusterizer::process_event(PHCompositeNode* topNode)
{
  // get node containing the digitized hits
  if (!do_read_raw)
  {
    m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    if (!m_hits)
    {
      std::cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  else
  {
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
    std::cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTER." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for cluster hit associations
  m_clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!m_clusterhitassoc)
  {
    std::cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTERHITASSOC" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for cluster-crossing associations
  m_clustercrossingassoc = findNode::getClass<TrkrClusterCrossingAssoc>(topNode, "TRKR_CLUSTERCROSSINGASSOC");
  if (!m_clustercrossingassoc)
  {
    std::cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTERCROSINGASSOC" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!do_read_raw)
  {
    ClusterLadderCells(topNode);
  }
  else
  {
    ClusterLadderCellsRaw(topNode);
  }
  PrintClusters(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

void InttClusterizer::CalculateLadderThresholds(PHCompositeNode* topNode)
{
  /*
  PHG4CellContainer* cells = findNode::getClass<PHG4CellContainer>(topNode, "G4CELL_INTT");
  if (!cells) return;
  */

  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if (!geom_container)
  {
    return;
  }

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
      _make_z_clustering.insert(std::make_pair(layer, false));
    }

    if (_make_e_weights.find(layer) == _make_e_weights.end())
    {
      _make_e_weights.insert(std::make_pair(layer, false));
    }
  }

  return;
}

void InttClusterizer::ClusterLadderCells(PHCompositeNode* topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "Entering InttClusterizer::ClusterLadderCells " << std::endl;
  }

  //----------
  // Get Nodes
  //----------

  // get the geometry node
  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if (!geom_container)
  {
    return;
  }

  //-----------
  // Clustering
  //-----------

  // loop over the InttHitSet objects
  TrkrHitSetContainer::ConstRange hitsetrange =
      m_hits->getHitSets(TrkrDefs::TrkrId::inttId);
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {
    // Each hitset contains only hits that are clusterizable - i.e. belong to a single sensor
    TrkrHitSet* hitset = hitsetitr->second;

    if (Verbosity() > 1)
    {
      std::cout << "InttClusterizer found hitsetkey " << hitsetitr->first << std::endl;
    }
    if (Verbosity() > 2)
    {
      hitset->identify();
    }

    // we have a single hitset, get the info that identifies the sensor
    int layer = TrkrDefs::getLayer(hitsetitr->first);
    int ladder_z_index = InttDefs::getLadderZId(hitsetitr->first);

    // we will need the geometry object for this layer to get the global position
    CylinderGeomIntt* geom = dynamic_cast<CylinderGeomIntt*>(geom_container->GetLayerGeom(layer));
    float pitch = geom->get_strip_y_spacing();
    float length = geom->get_strip_z_spacing();

    // fill a vector of hits to make things easier - gets every hit in the hitset
    std::vector<std::pair<TrkrDefs::hitkey, TrkrHit*>> hitvec;
    TrkrHitSet::ConstRange hitrangei = hitset->getHits();
    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
         hitr != hitrangei.second;
         ++hitr)
    {
      hitvec.emplace_back(hitr->first, hitr->second);
    }
    if (Verbosity() > 2)
    {
      std::cout << "hitvec.size(): " << hitvec.size() << std::endl;
    }

    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;
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
    std::vector<int> component(num_vertices(G));

    // this is the actual clustering, performed by boost
    connected_components(G, &component[0]);

    // Loop over the components(hit cells) compiling a list of the
    // unique connected groups (ie. clusters).
    std::set<int> cluster_ids;  // unique components

    std::multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>> clusters;
    for (unsigned int i = 0; i < component.size(); i++)
    {
      cluster_ids.insert(component[i]);                          // one entry per unique cluster id
      clusters.insert(std::make_pair(component[i], hitvec[i]));  // multiple entries per unique cluster id
    }

    // loop over the cluster ID's and make the clusters from the connected hits
    for (int clusid : cluster_ids)
    {
      // std::cout << " intt clustering: add cluster number " << clusid << std::endl;
      // get all hits for this cluster ID only
      std::pair<std::multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator,
                std::multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator>
          clusrange = clusters.equal_range(clusid);

      // make the cluster directly in the node tree
      TrkrDefs::cluskey ckey = TrkrDefs::genClusKey(hitset->getHitSetKey(), clusid);

      if (Verbosity() > 2)
      {
        std::cout << "Filling cluster with key " << ckey << std::endl;
      }

      // get the bunch crossing number from the hitsetkey
      short int crossing = InttDefs::getTimeBucketId(hitset->getHitSetKey());

      // determine the size of the cluster in phi and z, useful for track fitting the cluster
      std::set<int> phibins;
      std::set<int> zbins;

      // determine the cluster position...
      double xlocalsum = 0.0;
      double ylocalsum = 0.0;
      double zlocalsum = 0.0;
      unsigned int clus_adc = 0.0;
      unsigned int clus_maxadc = 0.0;
      unsigned nhits = 0;
      // std::cout << PHWHERE << " ckey " << ckey << ":" << std::endl;
      for (std::multimap<int, std::pair<TrkrDefs::hitkey, TrkrHit*>>::iterator mapiter = clusrange.first; mapiter != clusrange.second; ++mapiter)
      {
        // mapiter->second.first  is the hit key
        // std::cout << " adding hitkey " << mapiter->second.first << std::endl;
        int col = InttDefs::getCol((mapiter->second).first);
        int row = InttDefs::getRow((mapiter->second).first);
        zbins.insert(col);
        phibins.insert(row);

        // mapiter->second.second is the hit
        unsigned int hit_adc = (mapiter->second).second->getAdc();

        // Add clusterkey/bunch crossing to mmap
        m_clustercrossingassoc->addAssoc(ckey, crossing);

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
        if (hit_adc > clus_maxadc)
        {
          clus_maxadc = hit_adc;
        }
        clus_adc += hit_adc;
        ++nhits;

        // add this cluster-hit association to the association map of (clusterkey,hitkey)
        m_clusterhitassoc->addAssoc(ckey, mapiter->second.first);

        if (Verbosity() > 2)
        {
          std::cout << "     nhits = " << nhits << std::endl;
        }
        if (Verbosity() > 2)
        {
          std::cout << "  From  geometry object: hit x " << local_hit_location[0] << " hit y " << local_hit_location[1] << " hit z " << local_hit_location[2] << std::endl;
          std::cout << "     nhits " << nhits << " clusx  = " << xlocalsum / nhits << " clusy " << ylocalsum / nhits << " clusz " << zlocalsum / nhits << " hit_adc " << hit_adc << std::endl;
        }
      }

      static const float invsqrt12 = 1. / sqrt(12);

      // scale factors (phi direction)
      /*
        they corresponds to clusters of size 1 and 2 in phi
        other clusters, which are very few and pathological, get a scale factor of 1
        These scale factors are applied to produce cluster pulls with width unity
      */

      float phierror = pitch * invsqrt12;

      static constexpr std::array<double, 3> scalefactors_phi = {{0.85, 0.4, 0.33}};
      if (phibins.size() == 1 && layer < 5)
      {
        phierror *= scalefactors_phi[0];
      }
      else if (phibins.size() == 2 && layer < 5)
      {
        phierror *= scalefactors_phi[1];
      }
      else if (phibins.size() == 2 && layer > 4)
      {
        phierror *= scalefactors_phi[2];
      }
      // z error.
      const float zerror = zbins.size() * length * invsqrt12;

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

      auto clus = std::make_unique<TrkrClusterv5>();
      clus->setAdc(clus_adc);
      clus->setMaxAdc(clus_maxadc);
      clus->setLocalX(cluslocaly);
      clus->setLocalY(cluslocalz);
      clus->setPhiError(phierror);
      clus->setZError(zerror);
      clus->setPhiSize(phibins.size());
      clus->setZSize(zbins.size());
      // All silicon surfaces have a 1-1 map to hitsetkey.
      // So set subsurface key to 0
      clus->setSubSurfKey(0);

      if (Verbosity() > 2)
      {
        clus->identify();
      }

      m_clusterlist->addClusterSpecifyKey(ckey, clus.release());

    }  // end loop over cluster ID's
  }    // end loop over hitsets

  if (Verbosity() > 2)
  {
    // check that the associations were written correctly
    std::cout << "After InttClusterizer, cluster-hit associations are:" << std::endl;
    m_clusterhitassoc->identify();
  }

  if (Verbosity() > 0)
  {
    std::cout << " Cluster-crossing associations are:" << std::endl;
    m_clustercrossingassoc->identify();
  }

  return;
}
void InttClusterizer::ClusterLadderCellsRaw(PHCompositeNode* topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "Entering InttClusterizer::ClusterLadderCells " << std::endl;
  }

  //----------
  // Get Nodes
  //----------

  // get the geometry node
  PHG4CylinderGeomContainer* geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if (!geom_container)
  {
    return;
  }

  //-----------
  // Clustering
  //-----------

  // loop over the InttHitSet objects
  RawHitSetContainer::ConstRange hitsetrange =
      m_rawhits->getHitSets(TrkrDefs::TrkrId::inttId);
  for (RawHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {
    // Each hitset contains only hits that are clusterizable - i.e. belong to a single sensor
    RawHitSet* hitset = hitsetitr->second;

    if (Verbosity() > 1)
    {
      std::cout << "InttClusterizer found hitsetkey " << hitsetitr->first << std::endl;
    }
    if (Verbosity() > 2)
    {
      hitset->identify();
    }

    // we have a single hitset, get the info that identifies the sensor
    int layer = TrkrDefs::getLayer(hitsetitr->first);
    int ladder_z_index = InttDefs::getLadderZId(hitsetitr->first);

    // we will need the geometry object for this layer to get the global position
    CylinderGeomIntt* geom = dynamic_cast<CylinderGeomIntt*>(geom_container->GetLayerGeom(layer));
    float pitch = geom->get_strip_y_spacing();
    float length = geom->get_strip_z_spacing();

    // fill a vector of hits to make things easier - gets every hit in the hitset
    std::vector<RawHit*> hitvec;
    // int sector = InttDefs::getLadderPhiId(hitsetitr->first);
    // int side = InttDefs::getLadderZId(hitsetitr->first);

    RawHitSet::ConstRange hitrangei = hitset->getHits();
    for (RawHitSet::ConstIterator hitr = hitrangei.first;
         hitr != hitrangei.second;
         ++hitr)
    {
      // unsigned short iphi = (*hitr)->getPhiBin();
      // unsigned short it = (*hitr)->getTBin();
      //	std::cout << " intt layer " << layer << " sector: " << sector << " side " << side << " col: " << iphi << " row " << it << std::endl;
      hitvec.push_back((*hitr));
    }
    if (Verbosity() > 2)
    {
      std::cout << "hitvec.size(): " << hitvec.size() << std::endl;
    }

    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;
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
    std::vector<int> component(num_vertices(G));

    // this is the actual clustering, performed by boost
    connected_components(G, &component[0]);

    // Loop over the components(hit cells) compiling a list of the
    // unique connected groups (ie. clusters).
    std::set<int> cluster_ids;  // unique components

    std::multimap<int, RawHit*> clusters;
    for (unsigned int i = 0; i < component.size(); i++)
    {
      cluster_ids.insert(component[i]);                          // one entry per unique cluster id
      clusters.insert(std::make_pair(component[i], hitvec[i]));  // multiple entries per unique cluster id
    }

    // loop over the cluster ID's and make the clusters from the connected hits
    for (int clusid : cluster_ids)
    {
      // std::cout << " intt clustering: add cluster number " << clusid << std::endl;
      // get all hits for this cluster ID only
      auto clusrange = clusters.equal_range(clusid);

      // make the cluster directly in the node tree
      TrkrDefs::cluskey ckey = TrkrDefs::genClusKey(hitset->getHitSetKey(), clusid);

      if (Verbosity() > 2)
      {
        std::cout << "Filling cluster with key " << ckey << std::endl;
      }

      // get the bunch crossing number from the hitsetkey
      short int crossing = InttDefs::getTimeBucketId(hitset->getHitSetKey());

      // determine the size of the cluster in phi and z, useful for track fitting the cluster
      std::set<int> phibins;
      std::set<int> zbins;

      // determine the cluster position...
      double xlocalsum = 0.0;
      double ylocalsum = 0.0;
      double zlocalsum = 0.0;
      unsigned int clus_adc = 0.0;
      unsigned nhits = 0;

      // std::cout << PHWHERE << " ckey " << ckey << ":" << std::endl;

      std::map<int, unsigned int> m_phi, m_z;  // hold data for
      for (auto mapiter = clusrange.first; mapiter != clusrange.second; ++mapiter)
      {
        // mapiter->second.first  is the hit key
        // std::cout << " adding hitkey " << mapiter->second.first << std::endl;
        const auto energy = (mapiter->second)->getAdc();
        int col = (mapiter->second)->getPhiBin();
        int row = (mapiter->second)->getTBin();
        //	    std::cout << " found Tbin(row) " << row << " Phibin(col) " << col << std::endl;
        zbins.insert(col);
        phibins.insert(row);

        if (mClusHitsVerbose)
        {
          auto pnew = m_phi.try_emplace(row, energy);
          if (!pnew.second)
          {
            pnew.first->second += energy;
          }

          pnew = m_z.try_emplace(col, energy);
          if (!pnew.second)
          {
            pnew.first->second += energy;
          }
        }

        // mapiter->second.second is the hit
        unsigned int hit_adc = (mapiter->second)->getAdc();

        // Add clusterkey/bunch crossing to mmap
        m_clustercrossingassoc->addAssoc(ckey, crossing);

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
        //	    m_clusterhitassoc->addAssoc(ckey, mapiter->second.first);

        if (Verbosity() > 2)
        {
          std::cout << "     nhits = " << nhits << std::endl;
        }
        if (Verbosity() > 2)
        {
          std::cout << "  From  geometry object: hit x " << local_hit_location[0] << " hit y " << local_hit_location[1] << " hit z " << local_hit_location[2] << std::endl;
          std::cout << "     nhits " << nhits << " clusx  = " << xlocalsum / nhits << " clusy " << ylocalsum / nhits << " clusz " << zlocalsum / nhits << " hit_adc " << hit_adc << std::endl;
        }
      }

      if (mClusHitsVerbose)
      {
        if (Verbosity() > 10)
        {
          for (auto const& hit : m_phi)
          {
            std::cout << " m_phi(" << hit.first << " : " << hit.second << ") " << std::endl;
          }
        }
        for (auto& hit : m_phi)
        {
          mClusHitsVerbose->addPhiHit(hit.first, (float) hit.second);
        }
        for (auto& hit : m_z)
        {
          mClusHitsVerbose->addZHit(hit.first, (float) hit.second);
        }
        mClusHitsVerbose->push_hits(ckey);
      }

      static const float invsqrt12 = 1. / sqrt(12);

      // scale factors (phi direction)
      /*
        they corresponds to clusters of size 1 and 2 in phi
        other clusters, which are very few and pathological, get a scale factor of 1
        These scale factors are applied to produce cluster pulls with width unity
      */

      float phierror = pitch * invsqrt12;

      static constexpr std::array<double, 3> scalefactors_phi = {{0.85, 0.4, 0.33}};
      if (phibins.size() == 1 && layer < 5)
      {
        phierror *= scalefactors_phi[0];
      }
      else if (phibins.size() == 2 && layer < 5)
      {
        phierror *= scalefactors_phi[1];
      }
      else if (phibins.size() == 2 && layer > 4)
      {
        phierror *= scalefactors_phi[2];
      }
      // z error.
      const float zerror = zbins.size() * length * invsqrt12;

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

      auto clus = std::make_unique<TrkrClusterv5>();
      clus->setAdc(clus_adc);
      clus->setLocalX(cluslocaly);
      clus->setLocalY(cluslocalz);
      clus->setPhiError(phierror);
      clus->setZError(zerror);
      clus->setPhiSize(phibins.size());
      clus->setZSize(zbins.size());
      // All silicon surfaces have a 1-1 map to hitsetkey.
      // So set subsurface key to 0
      clus->setSubSurfKey(0);

      if (Verbosity() > 2)
      {
        clus->identify();
      }

      m_clusterlist->addClusterSpecifyKey(ckey, clus.release());

    }  // end loop over cluster ID's
  }    // end loop over hitsets

  if (Verbosity() > 2)
  {
    // check that the associations were written correctly
    std::cout << "After InttClusterizer, cluster-hit associations are:" << std::endl;
    m_clusterhitassoc->identify();
  }

  if (Verbosity() > 0)
  {
    std::cout << " Cluster-crossing associations are:" << std::endl;
    m_clustercrossingassoc->identify();
  }

  return;
}

void InttClusterizer::PrintClusters(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    TrkrClusterContainer* clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    if (!clusterlist)
    {
      return;
    }

    std::cout << "================= After InttClusterizer::process_event() ====================" << std::endl;

    std::cout << " There are " << clusterlist->size() << " clusters recorded: " << std::endl;

    clusterlist->identify();

    std::cout << "===========================================================================" << std::endl;
  }

  return;
}
