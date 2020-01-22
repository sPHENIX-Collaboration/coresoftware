/*!
 *  \file PHCASeeding.C
 *  \brief Progressive pattern recgnition based on GenFit Kalman filter
 *  \detail using Alan Dion's HelixHough for seeding, GenFit Kalman filter to do track propagation
 *  \author Christof Roland & Haiwang Yu
 */

//begin

#include "PHCASeeding.h"
#include "GPUTPCTrackLinearisation.h"
#include "GPUTPCTrackParam.h"

// trackbase_historic includes
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>  // for getLayer, clu...

// sPHENIX Geant4 includes
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

// sPHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHTimer.h>  // for PHTimer
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

//ROOT includes for debugging
#include <TFile.h>
#include <TNtuple.h>
#include <TVector3.h>  // for TVector3

//BOOST for combi seeding
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <utility>  // for pair, make_pair
#include <vector>
#include <algorithm> // for find

// forward declarations
class PHCompositeNode;



//#define _DEBUG_

#if defined(_DEBUG_)
#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#else
#define LogDebug(exp) (void)0
#endif

#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

//#define _DEBUG_

//end

typedef bg::model::point<float, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<point, TrkrDefs::cluskey> pointKey;
typedef std::array<TrkrDefs::cluskey,2> keylink;
typedef std::vector<TrkrDefs::cluskey> keylist;

using namespace std;
//using namespace ROOT::Minuit2;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

PHCASeeding::PHCASeeding(
    const string &name,
    unsigned int nlayers_maps,
    unsigned int nlayers_intt,
    unsigned int nlayers_tpc,
    unsigned int start_layer,
    float cluster_z_error,
    float cluster_alice_y_error,
    float neighbor_phi_width,
    float neighbor_eta_width,
    float maxSinPhi,
    float Bz)
  : PHTrackSeeding(name)
  , _g4clusters(nullptr)
  , _g4tracks(nullptr)
  , _g4vertexes(nullptr)
  , _cluster_map(nullptr)
  , _svtxhitsmap(nullptr)
  , _hit_used_map(nullptr)
  , _hit_used_map_size(0)
  , phisr(0.005)
  , etasr(0.0035)
  , phist(0.001)
  , etast(0.003)
  , phixt(0.008)
  , etaxt(0.005)
  , _vertex(nullptr)
  , _nlayers_maps(nlayers_maps)
  , _nlayers_intt(nlayers_intt)
  , _nlayers_tpc(nlayers_tpc)
  , _start_layer(start_layer)
  , _cluster_z_error(cluster_z_error)
  , _cluster_alice_y_error(cluster_alice_y_error)
  , _neighbor_phi_width(neighbor_phi_width)
  , _neighbor_eta_width(neighbor_eta_width)
  , _max_sin_phi(maxSinPhi)
  , _Bz(Bz)
  , _phi_scale(2)
  , _z_scale(2)
{
}

int PHCASeeding::InitializeGeometry(PHCompositeNode *topNode)
{
  PHG4CylinderCellGeomContainer *cellgeos = findNode::getClass<
      PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  PHG4CylinderGeomContainer *laddergeos = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  PHG4CylinderGeomContainer *mapsladdergeos = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");

  //_nlayers_seeding = _seeding_layer.size();
  //_radii.assign(_nlayers_seeding, 0.0);
  map<float, int> radius_layer_map;

  _radii_all.assign(60, 0.0);
  _layer_ilayer_map.clear();
  _layer_ilayer_map_all.clear();
  if (cellgeos)
  {
    PHG4CylinderCellGeomContainer::ConstRange layerrange =
        cellgeos->get_begin_end();
    for (PHG4CylinderCellGeomContainer::ConstIterator layeriter =
             layerrange.first;
         layeriter != layerrange.second; ++layeriter)
    {
      radius_layer_map.insert(
          make_pair(layeriter->second->get_radius(),
                    layeriter->second->get_layer()));
    }
  }

  if (laddergeos)
  {
    PHG4CylinderGeomContainer::ConstRange layerrange =
        laddergeos->get_begin_end();
    for (PHG4CylinderGeomContainer::ConstIterator layeriter =
             layerrange.first;
         layeriter != layerrange.second; ++layeriter)
    {
      radius_layer_map.insert(
          make_pair(layeriter->second->get_radius(),
                    layeriter->second->get_layer()));
    }
  }

  if (mapsladdergeos)
  {
    PHG4CylinderGeomContainer::ConstRange layerrange =
        mapsladdergeos->get_begin_end();
    for (PHG4CylinderGeomContainer::ConstIterator layeriter =
             layerrange.first;
         layeriter != layerrange.second; ++layeriter)
    {
      radius_layer_map.insert(
          make_pair(layeriter->second->get_radius(),
                    layeriter->second->get_layer()));
    }
  }
  for (map<float, int>::iterator iter = radius_layer_map.begin();
       iter != radius_layer_map.end(); ++iter)
  {
    _layer_ilayer_map_all.insert(make_pair(iter->second, _layer_ilayer_map_all.size()));

    /*if (std::find(_seeding_layer.begin(), _seeding_layer.end(),
                  iter->second) != _seeding_layer.end())
    {
      _layer_ilayer_map.insert(make_pair(iter->second, ilayer));
      ++ilayer;
      }*/
  }
  if (cellgeos)
  {
    PHG4CylinderCellGeomContainer::ConstRange begin_end =
        cellgeos->get_begin_end();
    PHG4CylinderCellGeomContainer::ConstIterator miter = begin_end.first;
    for (; miter != begin_end.second; ++miter)
    {
      PHG4CylinderCellGeom *geo = miter->second;
      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
          geo->get_radius() + 0.5 * geo->get_thickness();

      /*if (_layer_ilayer_map.find(geo->get_layer()) != _layer_ilayer_map.end())
      {
        _radii[_layer_ilayer_map[geo->get_layer()]] =
            geo->get_radius();
	    }*/
    }
  }

  if (laddergeos)
  {
    PHG4CylinderGeomContainer::ConstRange begin_end =
        laddergeos->get_begin_end();
    PHG4CylinderGeomContainer::ConstIterator miter = begin_end.first;
    for (; miter != begin_end.second; ++miter)
    {
      PHG4CylinderGeom *geo = miter->second;
      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
          geo->get_radius() + 0.5 * geo->get_thickness();

      /*if (_layer_ilayer_map.find(geo->get_layer()) != _layer_ilayer_map.end())
      {
        _radii[_layer_ilayer_map[geo->get_layer()]] = geo->get_radius();
	}*/
    }
  }

  if (mapsladdergeos)
  {
    PHG4CylinderGeomContainer::ConstRange begin_end =
        mapsladdergeos->get_begin_end();
    PHG4CylinderGeomContainer::ConstIterator miter = begin_end.first;
    for (; miter != begin_end.second; ++miter)
    {
      PHG4CylinderGeom *geo = miter->second;

      //if(geo->get_layer() > (int) _radii.size() ) continue;

      //			if (Verbosity() >= 2)
      //				geo->identify();

      //TODO
      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
          geo->get_radius();

      /*if (_layer_ilayer_map.find(geo->get_layer()) != _layer_ilayer_map.end())
      {
        _radii[_layer_ilayer_map[geo->get_layer()]] = geo->get_radius();
	}*/
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}


int PHCASeeding::GetNodes(PHCompositeNode *topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

double PHCASeeding::phiadd(double phi1, double phi2)
{
  double s = phi1 + phi2;
  if (s > 2 * M_PI)
    return s - 2 * M_PI;
  else if (s < 0)
    return s + 2 * M_PI;
  else
    return s;
}

double PHCASeeding::phidiff(double phi1, double phi2)
{
  double d = phi1 - phi2;
  if (d > M_PI)
    return d - 2 * M_PI;
  else if (d < -M_PI)
    return d + 2 * M_PI;
  else
    return d;
}

void PHCASeeding::QueryTree(const bgi::rtree<pointKey, bgi::quadratic<16>> &rtree, double phimin, double etamin, double lmin, double phimax, double etamax, double lmax, std::vector<pointKey> &returned_values)
{
  rtree.query(bgi::intersects(box(point(phimin, etamin, lmin), point(phimax, etamax, lmax))), std::back_inserter(returned_values));
  if (phimin < 0) rtree.query(bgi::intersects(box(point(2 * M_PI + phimin, etamin, lmin), point(2 * M_PI, etamax, lmax))), std::back_inserter(returned_values));
  if (phimax > 2 * M_PI) rtree.query(bgi::intersects(box(point(0, etamin, lmin), point(phimax - 2 * M_PI, etamax, lmax))), std::back_inserter(returned_values));
}

void PHCASeeding::FillTree()
{
  //  bgi::rtree<pointKey, bgi::quadratic<16> > rtree;
  PHTimer *t_fill = new PHTimer("t_fill");
  t_fill->stop();
  int n_dupli = 0;
  int nlayer[60];
  for (int j = 0; j < 60; ++j) nlayer[j] = 0;

  TrkrClusterContainer::ConstRange clusrange = _cluster_map->getClusters();

  for (TrkrClusterContainer::ConstIterator iter = clusrange.first; iter != clusrange.second; ++iter)
  {
    TrkrCluster *cluster = iter->second;
    TrkrDefs::cluskey ckey = iter->first;
    unsigned int layer = TrkrDefs::getLayer(ckey);
    if (layer < (_nlayers_maps + _nlayers_intt)) continue;

    TVector3 vec(cluster->getPosition(0) - _vertex->get_x(), cluster->getPosition(1) - _vertex->get_y(), cluster->getPosition(2) - _vertex->get_z());

    double clus_phi = vec.Phi();
    clus_phi -= 2 * M_PI * floor(clus_phi / (2 * M_PI));
    double clus_eta = vec.Eta();
    double clus_l = layer;  // _radii_all[layer];

    vector<pointKey> testduplicate;
    QueryTree(_rtree, clus_phi - 0.00001, clus_eta - 0.00001, layer - 0.5, clus_phi + 0.00001, clus_eta + 0.00001, layer + 0.5, testduplicate);
    if (!testduplicate.empty())
    {
      ++n_dupli;
      continue;
    }
    ++nlayer[layer];
    t_fill->restart();
    _rtree.insert(std::make_pair(point(clus_phi, clus_eta, clus_l), ckey));
    t_fill->stop();
  }

  std::cout << "fill time: " << t_fill->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "number of duplicates : " << n_dupli << std::endl;
}

int PHCASeeding::Process(PHCompositeNode *topNode)
{
  TFile fpara("CA_para.root", "RECREATE");
  TNtuple *NT = new TNtuple("NT", "NT", "pt:dpt:z:dz:phi:dphi:c:dc:nhit");
  _vertex = _vertex_map->get(0);

  //for different purpose
  phisr = 0.005 * _phi_scale;  // *2
  etasr = 0.0035 * _z_scale;   // *2;
  /* 0.7 version
     0.9 version 
     phist = 0.001*2;
     etast = 0.003*2;
  */
  phist = 0.001 * _phi_scale;  // *5 *7;
  etast = 0.003 * _z_scale;    // *5/ *7;

  PHTimer *t_seed = new PHTimer("t_seed");
  t_seed->stop();
  t_seed->restart();

  _rtree.clear();
  FillTree();

  int numberofseeds = 0;
  LogDebug(" entries in tree: " << _rtree.size() << endl);

  vector<pointKey> allClusters;
  vector<keylink> belowLinks;
  vector<keylink> aboveLinks;
  // messy way of getting vector<pointKey> for all clusters
  QueryTree(_rtree,
            0, // phi
            -3, // eta
            -1, // layer 
            2*M_PI, // phi
            3, // eta
            100, // layer
            allClusters);
  LogDebug(" number of total clusters: " << allClusters.size() << endl);
  for (vector<pointKey>::iterator StartCluster = allClusters.begin(); StartCluster != allClusters.end(); ++StartCluster)
  {
    // get clusters near this one in adjacent layers
    double StartPhi = StartCluster->first.get<0>();
    double StartEta = StartCluster->first.get<1>();
    unsigned int StartLayer = TrkrDefs::getLayer(StartCluster->second);
    //if(StartLayer > _start_layer) continue;
    TrkrCluster* StartCl = _cluster_map->findCluster(StartCluster->second);
    double StartX = StartCl->getPosition(0);
    double StartY = StartCl->getPosition(1);
    double StartZ = StartCl->getPosition(2);
    LogDebug(" starting cluster:" << endl);
    LogDebug(" eta: " << StartEta << endl);
    LogDebug(" phi: " << StartPhi << endl);
    LogDebug(" layer: " << StartLayer << endl);

    vector<pointKey> ClustersAbove;
    vector<pointKey> ClustersBelow;
    QueryTree(_rtree,
              StartPhi-_neighbor_phi_width,
              StartEta-_neighbor_eta_width,
              (double) StartLayer - 1.5,
              StartPhi+_neighbor_phi_width,
              StartEta+_neighbor_eta_width,
              (double) StartLayer - 0.5,
              ClustersBelow);
    QueryTree(_rtree,
              StartPhi-_neighbor_phi_width,
              StartEta-_neighbor_eta_width,
              (double) StartLayer + 0.5,
              StartPhi+_neighbor_phi_width,
              StartEta+_neighbor_eta_width,
              (double) StartLayer + 1.5,
              ClustersAbove);
    LogDebug(" entries in below layer: " << ClustersBelow.size() << endl);
    LogDebug(" entries in above layer: " << ClustersAbove.size() << endl);
    vector<array<double,3>> delta_below;
    vector<array<double,3>> delta_above;
    // calculate (delta_eta, delta_phi) vector for each neighboring cluster
    for (vector<pointKey>::iterator BelowCandidate = ClustersBelow.begin(); BelowCandidate != ClustersBelow.end(); ++BelowCandidate)
    {
      TrkrCluster* BelowCl = _cluster_map->findCluster(BelowCandidate->second);
      double BelowX = BelowCl->getPosition(0)-StartX;
      double BelowY = BelowCl->getPosition(1)-StartY;
      double BelowZ = BelowCl->getPosition(2)-StartZ;
      delta_below.push_back(
        {BelowX,
         BelowY,
         BelowZ});
    }
    for(vector<pointKey>::iterator AboveCandidate = ClustersAbove.begin(); AboveCandidate != ClustersAbove.end(); ++AboveCandidate)
    {
      TrkrCluster* AboveCl = _cluster_map->findCluster(AboveCandidate->second);
      double AboveX = AboveCl->getPosition(0)-StartX;
      double AboveY = AboveCl->getPosition(1)-StartY;
      double AboveZ = AboveCl->getPosition(2)-StartZ;
      delta_above.push_back(
        {AboveX,
         AboveY,
         AboveZ});
    }
    // find the three clusters closest to a straight line
    // (by maximizing the cos of the angle between the (delta_eta,delta_phi) vectors)
    double maxCosPlaneAngle = 0;
    TrkrDefs::cluskey bestBelowCluster = 0;
    TrkrDefs::cluskey bestAboveCluster = 0;
    for(size_t iAbove = 0; iAbove<delta_above.size(); ++iAbove)
    {
      for(size_t iBelow = 0; iBelow<delta_below.size(); ++iBelow)
      {
        double dotProduct = delta_below[iBelow][0]*delta_above[iAbove][0]+delta_below[iBelow][1]*delta_above[iAbove][1]+delta_below[iBelow][2]*delta_above[iAbove][2];
        double belowSqLength = sqrt(pow(delta_below[iBelow][0],2)+pow(delta_below[iBelow][1],2)+pow(delta_below[iBelow][2],2));
        double aboveSqLength = sqrt(pow(delta_above[iAbove][0],2)+pow(delta_above[iAbove][1],2)+pow(delta_above[iAbove][2],2));
        double cosPlaneAngle = dotProduct / (belowSqLength*aboveSqLength);
        if(cosPlaneAngle < maxCosPlaneAngle )
        {
          maxCosPlaneAngle = cosPlaneAngle;
          bestBelowCluster = ClustersBelow[iBelow].second;
          bestAboveCluster = ClustersAbove[iAbove].second;
        }
      }
    }
    belowLinks.push_back({StartCluster->second,bestBelowCluster});
    aboveLinks.push_back({StartCluster->second,bestAboveCluster});
    LogDebug(" max collinearity: " << maxCosPlaneAngle << endl);
#if defined(_DEBUG_)
    if(bestBelowCluster==0 || bestAboveCluster == 0)
    {
      LogDebug("Incomplete triplet, skipping debug output" << endl);
    }
    else
    {
      TrkrCluster* bcl = _cluster_map->findCluster(bestBelowCluster);
      TrkrCluster* scl = _cluster_map->findCluster(StartCluster->second);
      TrkrCluster* acl = _cluster_map->findCluster(bestAboveCluster);
      LogDebug(" found triplet: (" << bcl->getPosition(0) << "," << bcl->getPosition(1) << "," << bcl->getPosition(2)
       << ")<-(" << scl->getPosition(0) << "," << scl->getPosition(1) << "," << scl->getPosition(2) << ")->(" 
       << acl->getPosition(0) << "," << acl->getPosition(1) << "," << acl->getPosition(2) << ")" << endl);
    }
#endif
  }
  // remove all triplets for which there isn't a mutual association between two clusters
  vector<keylink> bidirectionalLinks;
  for(vector<keylink>::iterator belowLink = belowLinks.begin(); belowLink != belowLinks.end(); ++belowLink)
  {
    keylink reversed = {(*belowLink)[1],(*belowLink)[0]};
    vector<keylink>::iterator sameAboveLinkExists = find(aboveLinks.begin(),aboveLinks.end(),reversed);
    if(sameAboveLinkExists != aboveLinks.end())
    {
      bidirectionalLinks.push_back((*belowLink));
    }
  }
  LogDebug(" bidirectional links found:" << endl);
  for(vector<keylink>::iterator l = bidirectionalLinks.begin(); l != bidirectionalLinks.end(); ++l)
  {
    LogDebug(" "<< (*l)[0] << " <-> " << (*l)[1] << endl);
  }
  // follow bidirectional links to form lists of cluster keys
  // (to be fitted for track seed parameters)
  vector<keylist> trackSeedKeyLists;
  // get starting cluster keys, create a keylist for each
  // (only check last element of each pair because we start from the outer layers and go inward)
  for(vector<keylink>::iterator startCand = bidirectionalLinks.begin(); startCand != bidirectionalLinks.end(); ++startCand)
  {
    bool has_above_link = false;
    for(vector<keylink>::iterator testlink = bidirectionalLinks.begin(); testlink != bidirectionalLinks.end(); ++testlink)
    {
      if((*startCand) == (*testlink)) continue;
      if((*startCand)[0] == (*testlink)[1]) has_above_link = true;
    }
    if(!has_above_link)
    {
      trackSeedKeyLists.push_back({(*startCand)[0],(*startCand)[1]});
    }
  }
  // assemble track cluster chains from starting cluster keys (ordered from outside in)
  for(vector<keylist>::iterator trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
  {
    bool reached_end = false;
    while(!reached_end)
    {
      TrkrDefs::cluskey trackHead = trackKeyChain->back();
      bool no_next_link = true;
      for(vector<keylink>::iterator testlink = bidirectionalLinks.begin(); testlink != bidirectionalLinks.end(); ++testlink)
      {
        if((*testlink)[0]==trackHead)
        {
          trackKeyChain->push_back((*testlink)[1]);
          no_next_link = false;
        }
      }
      if(no_next_link) reached_end = true;
    }
  }
  LogDebug(" track key chains assembled: " << trackSeedKeyLists.size() << endl);
  LogDebug(" track key chain lengths: " << endl);
  for(vector<keylist>::iterator trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
  {
    LogDebug(" " << trackKeyChain->size() << endl);
  }
  int jumpcount = 0;
  LogDebug(" track key associations:" << endl);
  for(size_t i=0;i<trackSeedKeyLists.size();++i)
  {
    LogDebug(" seed " << i << ":" << endl);
    double lasteta = -100;
    double lastphi = -100;
    for(size_t j=0;j<trackSeedKeyLists[i].size();++j)
    {
      TrkrCluster* cl = _cluster_map->findCluster(trackSeedKeyLists[i][j]);
      TVector3 vec(cl->getPosition(0) - _vertex->get_x(), cl->getPosition(1) - _vertex->get_y(), cl->getPosition(2) - _vertex->get_z());

      double clus_phi = vec.Phi();
      clus_phi -= 2 * M_PI * floor(clus_phi / (2 * M_PI));
      double clus_eta = vec.Eta();
      double etajump = clus_eta-lasteta;
      double phijump = clus_phi-lastphi;
      #if defined(_DEBUG_) 
      unsigned int lay = TrkrDefs::getLayer(trackSeedKeyLists[i][j]);
      #endif
      if((fabs(etajump)>0.1 && lasteta!=-100) || (fabs(phijump)>1 && lastphi!=-100))
      {
         LogDebug(" Eta or Phi jump too large! " << endl);
         ++jumpcount;
      }
      LogDebug(" (eta,phi,layer) = (" << clus_eta << "," << clus_phi << "," << lay << ") " <<
        " (x,y,z) = (" << cl->getPosition(0) << "," << cl->getPosition(1) << "," << cl->getPosition(2) << ")" << endl);
      lasteta = clus_eta;
      lastphi = clus_phi;
    }
  }
  LogDebug(" Total large jumps: " << jumpcount << endl);
  // Turn track cluster chains into track candidates using ALICE simplified KF.

  for(vector<keylist>::iterator trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
  {
    // get starting cluster from key
    TrkrCluster* startCluster = _cluster_map->findCluster(trackKeyChain->at(0));
    // Transform sPHENIX coordinates into ALICE-compatible coordinates
    double x0 = startCluster->getPosition(0);
    double y0 = startCluster->getPosition(1);
    double z0 = startCluster->getPosition(2);
    LogDebug("Initial (x,y,z): (" << x0 << "," << y0 << "," << z0 << ")" << endl);
    // ALICE x coordinate = distance from beampipe
    float alice_x0 = sqrt(x0*x0+y0*y0);
    float alice_y0 = 0;
    float alice_z0 = z0;
    // Initialize track and linearisation
    GPUTPCTrackParam trackSeed;
    trackSeed.InitParam();
    trackSeed.SetX(alice_x0);
    trackSeed.SetY(alice_y0);
    trackSeed.SetZ(alice_z0);
    float x = x0;
    float y = y0;
    #if defined(_DEBUG_)
    float z = z0;
    float alice_x = sqrt(x0*x0+y0*y0);
    #endif
    float trackCartesian_x = 0.;
    float trackCartesian_y = 0.;
    float trackCartesian_z = 0.;
    // Pre-set momentum-based parameters to improve numerical stability
    TrkrCluster* SecondCluster = _cluster_map->findCluster(trackKeyChain->at(1));
    float second_x = SecondCluster->getPosition(0);
    float second_y = SecondCluster->getPosition(1);
    float second_z = SecondCluster->getPosition(2);
    float second_alice_x = sqrt(second_x*second_x+second_y*second_y);
    float delta_alice_x = second_alice_x - alice_x0;
    float first_phi = atan(y0/x0);
    float second_alice_y = (second_x/cos(first_phi)-second_y/sin(first_phi))/(sin(first_phi)/cos(first_phi)+cos(first_phi)/sin(first_phi));
    float init_SinPhi = second_alice_y / sqrt(delta_alice_x*delta_alice_x + second_alice_y*second_alice_y);
    float delta_z = second_z - z0;
    float init_DzDs = delta_z / sqrt(delta_alice_x*delta_alice_x + second_alice_y*second_alice_y);
    trackSeed.SetSinPhi(init_SinPhi);
    LogDebug("Set initial SinPhi to " << init_SinPhi << endl);
    trackSeed.SetDzDs(init_DzDs);
    LogDebug("Set initial DzDs to " << init_DzDs << endl);
    GPUTPCTrackLinearisation trackLine(trackSeed);

    LogDebug(endl << endl << "------------------------" << endl << "seed size: " << trackKeyChain->size() << endl << endl << endl);
    int cluster_ctr = 1;
    // starting at second cluster, perform track propagation
    for(keylist::iterator clusterkey = next(trackKeyChain->begin()); clusterkey != trackKeyChain->end(); ++clusterkey)
    {
      LogDebug("cluster " << cluster_ctr << " -> " << cluster_ctr + 1 << endl);
      LogDebug("this cluster (x,y,z) = (" << x << "," << y << "," << z << ")" << endl);
      // get cluster from key
      TrkrCluster* nextCluster = _cluster_map->findCluster(*clusterkey);
      // find ALICE x-coordinate
      float nextCluster_x = nextCluster->getPosition(0);
      float nextCluster_y = nextCluster->getPosition(1);
      float nextCluster_z = nextCluster->getPosition(2);
      float nextAlice_x = sqrt(nextCluster_x*nextCluster_x+nextCluster_y*nextCluster_y);
      // rotate track coordinates to match orientation of next cluster
      float newPhi = atan(nextCluster_y/nextCluster_x);
      LogDebug("new phi = " << newPhi << endl);
      float oldPhi = atan(y/x);
      LogDebug("old phi = " << oldPhi << endl);
      float alpha = newPhi - oldPhi;
      LogDebug("alpha = " << alpha << endl);
      if(!trackSeed.Rotate(alpha,trackLine,_max_sin_phi))
      {
        LogError("Rotate failed! Aborting for this seed...");
        break;
      }
      LogDebug("track coordinates (ALICE) after rotation: (" << trackSeed.GetX() << "," << trackSeed.GetY() << "," << trackSeed.GetZ() << ")" << endl);
      LogDebug("Transporting from " << alice_x << " to " << nextAlice_x << "...");
      if(!trackSeed.TransportToX(nextAlice_x,trackLine,_Bz,_max_sin_phi))
      {
        LogError("Transport failed! Aborting for this seed...");
        break;
      }
      // convert ALICE coordinates to sPHENIX cartesian coordinates, for debugging
      float predicted_alice_x = trackSeed.GetX();
      LogDebug("new track ALICE x = " << trackSeed.GetX() << endl);
      float predicted_alice_y = trackSeed.GetY();
      LogDebug("new track ALICE y = " << trackSeed.GetY() << endl);
      float predicted_z = trackSeed.GetZ();
      LogDebug("new track z = " << trackSeed.GetZ() << endl);
      float cos_phi = x/sqrt(x*x+y*y);
      LogDebug("cos_phi = " << cos_phi << endl);
      float sin_phi = y/sqrt(x*x+y*y);
      LogDebug("sin phi = " << sin_phi << endl);
      trackCartesian_x = predicted_alice_x*cos_phi+predicted_alice_y*sin_phi;
      trackCartesian_y = predicted_alice_x*sin_phi-predicted_alice_y*cos_phi;
      trackCartesian_z = predicted_z;
      LogDebug("Track transported to (x,y,z) = (" << trackCartesian_x << "," << trackCartesian_y << "," << trackCartesian_z << ")" << endl);
      LogDebug("Next cluster is at (x,y,z) = (" << nextCluster_x << "," << nextCluster_y << "," << nextCluster_z << ")" << endl);
      LogDebug("track coordinates (ALICE) after rotation: (" << trackSeed.GetX() << "," << trackSeed.GetY() << "," << trackSeed.GetZ() << ")" << endl);
      //float nextCluster_alice_y = (nextCluster_x/cos(newPhi) - nextCluster_y/sin(newPhi))/(tan(newPhi)+1./tan(newPhi));
      float nextCluster_alice_y = 0.;
      LogDebug("next cluster ALICE y = " << nextCluster_alice_y << endl);
      float y2_error = _cluster_alice_y_error*_cluster_alice_y_error;
      float z2_error = _cluster_z_error*_cluster_z_error;
      LogDebug("track ALICE SinPhi = " << trackSeed.GetSinPhi() << endl);
      // Apply Kalman filter
      if(!trackSeed.Filter(nextCluster_alice_y,nextCluster_z,y2_error,z2_error,_max_sin_phi))
      {
        LogError("Kalman filter failed for seed " << numberofseeds << "! Aborting for this seed..." << endl);
        break;
      }
      x = nextCluster_x;
      y = nextCluster_y;
      #if defined(_DEBUG_)
      z = nextCluster_z;
      alice_x = nextAlice_x;
      #endif
      ++cluster_ctr;
    } 
    //    pt:z:dz:phi:dphi:c:dc
    // Fill NT with track parameters
    float StartEta = -log(tan(atan(z0/sqrt(x0*x0+y0*y0))));
    float track_pt = abs( 1./(trackSeed.GetQPt()));
    LogDebug("Track pt = " << track_pt << endl);
    float track_pterr = sqrt(trackSeed.GetErr2QPt())/(trackSeed.GetQPt()*trackSeed.GetQPt());
    LogDebug("Track pterr = " << track_pterr << endl);
    float track_z = trackSeed.GetZ();
    float track_zerr = sqrt(trackSeed.GetErr2Z());
    float track_phi = atan(trackCartesian_y/trackCartesian_x);
    float last_cluster_phierr = _cluster_map->findCluster(trackKeyChain->back())->getPhiError();
    // phi error assuming error in track radial coordinate is zero
    float track_phierr = sqrt(pow(last_cluster_phierr,2)+(pow(trackSeed.GetX(),2)*trackSeed.GetErr2Y()) / 
      pow(pow(trackSeed.GetX(),2)+pow(trackSeed.GetY(),2),2));
    LogDebug("Track phi = " << track_phi << endl);
    LogDebug("Track phierr = " << track_phierr << endl);
    float track_curvature = trackSeed.GetKappa(_Bz);
    float track_curverr = sqrt(trackSeed.GetErr2QPt())*_Bz;
    NT->Fill(track_pt, track_pterr, track_z, track_zerr, track_phi, track_phierr, track_curvature, track_curverr, trackKeyChain->size());
    SvtxTrack_v1 track;
    track.set_id(numberofseeds);
    for (unsigned int j = 0; j < trackKeyChain->size(); ++j)
    {
      track.insert_cluster_key(trackKeyChain->at(j));
    }
    track.set_ndf(trackSeed.GetNDF());
    if(trackSeed.GetQPt()<0) track.set_charge(-1);
    else track.set_charge(1);
    TrkrCluster *cl = _cluster_map->findCluster(trackKeyChain->at(0));
    track.set_x(trackCartesian_x);  //track.set_x(cl->getX());
    track.set_y(trackCartesian_y);  //track.set_y(cl->getY());
    track.set_z(trackCartesian_z);  //track.set_z(cl->getZ());
    track.set_px(track_pt * cos(track_phi));
    track.set_py(track_pt * sin(track_phi));
    track.set_pz(track_pt / tan(2 * atan(exp(-StartEta))));
    track.set_error(0, 0, cl->getError(0, 0));
    track.set_error(0, 1, cl->getError(0, 1));
    track.set_error(0, 2, cl->getError(0, 2));
    track.set_error(1, 1, cl->getError(1, 1));
    track.set_error(1, 2, cl->getError(1, 2));
    track.set_error(2, 2, cl->getError(2, 2));
    track.set_error(3, 3, track_pterr * track_pterr * cos(track_phi) * cos(track_phi));
    track.set_error(4, 4, track_pterr * track_pterr * sin(track_phi) * sin(track_phi));
    track.set_error(5, 5, track_pterr * track_pterr / tan(2 * atan(exp(-StartEta))) / tan(2 * atan(exp(-StartEta))));
    _track_map->insert(&track);
    ++numberofseeds;
  }
  t_seed->stop();
  cout << "number of seeds " << numberofseeds << endl;
  cout << "seeding time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  fpara.cd();
  NT->Write();
  fpara.Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCASeeding::Setup(PHCompositeNode *topNode)
{
  cout << "Called Setup" << endl;
  cout << "topNode:" << topNode << endl;
  PHTrackSeeding::Setup(topNode);
  //  int ret = GetNodes(topNode);
  //return ret;
  GetNodes(topNode);
  InitializeGeometry(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCASeeding::End()
{
  cout << "Called End " << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
