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

//end

typedef bg::model::point<float, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<point, TrkrDefs::cluskey> pointKey;
typedef std::array<TrkrDefs::cluskey,2> keylink;
typedef std::vector<TrkrDefs::cluskey> keylist;

#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

using namespace std;
//using namespace ROOT::Minuit2;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

PHCASeeding::PHCASeeding(
    const string &name,
    unsigned int nlayers_maps,
    unsigned int nlayers_intt,
    unsigned int nlayers_tpc,
    unsigned int start_layer)
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
  , _nlayers_maps(nlayers_maps)
  , _nlayers_intt(nlayers_intt)
  , _nlayers_tpc(nlayers_tpc)
  , _start_layer(start_layer)
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
    for (; miter != begin_end.second; miter++)
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
    for (; miter != begin_end.second; miter++)
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
    for (; miter != begin_end.second; miter++)
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
  for (int j = 0; j < 60; j++) nlayer[j] = 0;

  TrkrClusterContainer::ConstRange clusrange = _cluster_map->getClusters();

  for (TrkrClusterContainer::ConstIterator iter = clusrange.first; iter != clusrange.second; ++iter)
  {
    TrkrCluster *cluster = iter->second;
    TrkrDefs::cluskey ckey = iter->first;
    unsigned int layer = TrkrDefs::getLayer(ckey);
    if (layer < 39) continue;

    TVector3 vec(cluster->getPosition(0) - _vertex->get_x(), cluster->getPosition(1) - _vertex->get_y(), cluster->getPosition(2) - _vertex->get_z());

    double clus_phi = vec.Phi();
    clus_phi -= 2 * M_PI * floor(clus_phi / (2 * M_PI));
    double clus_eta = vec.Eta();
    double clus_l = layer;  // _radii_all[layer];

    vector<pointKey> testduplicate;
    QueryTree(_rtree, clus_phi - 0.00001, clus_eta - 0.00001, layer - 0.5, clus_phi + 0.00001, clus_eta + 0.00001, layer + 0.5, testduplicate);
    if (!testduplicate.empty())
    {
      n_dupli++;
      continue;
    }
    nlayer[layer]++;
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
  for (vector<pointKey>::iterator StartCluster = allClusters.begin(); StartCluster != allClusters.end(); StartCluster++)
  {
    // get clusters near this one in adjacent layers
    double StartPhi = StartCluster->first.get<0>();
    double StartEta = StartCluster->first.get<1>();
    unsigned int StartLayer = TrkrDefs::getLayer(StartCluster->second);
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
              StartPhi-M_PI/6,
              StartEta-1,
              (double) StartLayer - 1.5,
              StartPhi+M_PI/6,
              StartEta+1,
              (double) StartLayer - 0.5,
              ClustersBelow);
    QueryTree(_rtree,
              StartPhi-M_PI/6,
              StartEta-1,
              (double) StartLayer + 0.5,
              StartPhi+M_PI/6,
              StartEta+1,
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
    cout << "am I crashing here?" << endl;
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
      unsigned int lay = TrkrDefs::getLayer(trackSeedKeyLists[i][j]);
      if((fabs(etajump)>0.1 && lasteta!=-100) || (fabs(phijump)>1 && lastphi!=-100))
      {
         LogDebug(" Eta or Phi jump too large! " << endl);
         jumpcount++;
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
    GPUTPCTrackLinearisation trackLine; // default constructor is fine for now
    // convert B field to ALICE-compatible units
    float Bz = 1.4*0.000299792458f;
    // configurable max sin phi (algorithm doesn't work when track is too horizontal)
    float maxSinPhi = 0.99;
    float x = x0;
    float y = y0;
    float z = z0;
    // starting at SECOND cluster, perform track propagation
    for(keylist::iterator clusterkey = next(trackKeyChain->begin()); clusterkey != trackKeyChain->end(); ++clusterkey)
    {
      // get cluster from key
      TrkrCluster* nextCluster = _cluster_map->findCluster(*clusterkey);
      // find ALICE x-coordinate, and transport to that x
      float nextCluster_x = nextCluster->getPosition(0);
      float nextCluster_y = nextCluster->getPosition(1);
      float nextCluster_z = nextCluster->getPosition(2);
      float nextAlice_x = sqrt(nextCluster_x*nextCluster_x+nextCluster_y*nextCluster_y);
      LogDebug("Transporting from " << x << " to " << nextCluster_x << "...");
      
      if(!trackSeed.TransportToX(nextAlice_x,trackLine,Bz,maxSinPhi))
      {
        LogError("Transport failed! Aborting for this seed...");
        break;
      }
      // convert ALICE coordinates to sPHENIX cartesian coordinates, for debugging
      float predicted_alice_y = trackSeed.GetY();
      float predicted_z = trackSeed.GetZ();
      float cos_theta = x/sqrt(x*x+y*y);
      float sin_theta = y/sqrt(x*x+y*y);
      float delta_y = predicted_alice_y*cos_theta;
      float delta_x = predicted_alice_y*sin_theta;
      float trackCartesian_x = x + delta_x;
      float trackCartesian_y = y + delta_y;
      float trackCartesian_z = predicted_z;
      cout << "Track transported to (x,y,z) = (" << trackCartesian_x << "," << trackCartesian_y << "," << trackCartesian_z << ")" << endl;
      cout << "Next cluster is at (x,y,z) = (" << nextCluster_x << "," << nextCluster_y << "," << nextCluster_z << ")" << endl;
      // Rotate track coordinate system to be parallel to layer
      float newTheta = atan(trackCartesian_y/trackCartesian_x);
      float oldTheta = atan(y/x);
      float alpha = newTheta - oldTheta;
      if(!trackSeed.Rotate(alpha,trackLine,maxSinPhi))
      {
        LogError("Rotate failed! Aborting for this seed...");
        break;
      }
      // Calculate squared errors in Y and Z
      float z2_error = (nextCluster_z-predicted_z)*(nextCluster_z-predicted_z);
      float nextCluster_alice_y = (nextCluster_x*x+nextCluster_y*y+nextCluster_z*z)/sqrt(x*x+y*y+z*z);
      float y2_error = (nextCluster_alice_y-predicted_alice_y)*(nextCluster_alice_y-predicted_alice_y);
      // Apply Kalman filter
      if(!trackSeed.Filter(nextCluster_alice_y,nextCluster_z,y2_error,z2_error,maxSinPhi))
      {
        LogError("Kalman filter failed! Aborting for this seed...");
        break;
      }
      x = trackCartesian_x;
      y = trackCartesian_y;
      z = trackCartesian_z;
    } 
/*
        for (unsigned int newlayer = _start_layer - 2; newlayer >= (_start_layer - 7); newlayer--)
        {
          vector<pointKey> newlayer_clusters;
          cout << " window - "
               << " phimin " << currentphi - dphidr * (_radii_all[lastgoodlayer] - _radii_all[newlayer]) - phist
               << " phimax " << currentphi - dphidr * (_radii_all[lastgoodlayer] - _radii_all[newlayer]) + phist
               << " etamin " << currenteta - etast
               << " etamax " << currenteta + etast
               << endl;
          QueryTree(_rtree,
                    currentphi - dphidr * (_radii_all[lastgoodlayer] - _radii_all[newlayer]) - phist,
                    currenteta - etast,
                    newlayer - 0.5,
                    currentphi - dphidr * (_radii_all[lastgoodlayer] - _radii_all[newlayer]) + phist,
                    currenteta + etast,
                    newlayer + 0.5,
                    newlayer_clusters);

          if (newlayer_clusters.empty())
          {
            failures += 1;
            if (failures > 2) break;  //0.7 2 0.9 3
          }
          else
          {
            double xinrecord = 100.0;
            pointKey *xinkey = &*newlayer_clusters.begin();
            for (std::vector<pointKey>::iterator it = newlayer_clusters.begin(); it != newlayer_clusters.end(); ++it)
            {
              double dist = abs(phidiff(it->first.get<0>(), currentphi - dphidr * (_radii_all[lastgoodlayer] - _radii_all[newlayer]))) + abs(it->first.get<1>() - currenteta);

              cout << " nuphi: " << it->first.get<0>()
                   << " nueta: " << it->first.get<1>()
                   << " dist: " << dist
                   << " lay: " << newlayer
                   << " dl: " << lastgoodlayer - newlayer
                   << " r: " << _radii_all[newlayer]
                   << " dr: " << _radii_all[lastgoodlayer] - _radii_all[newlayer]
                   << endl;
              if (dist < xinrecord)
              {
                *xinkey = *it;
                xinrecord = dist;
              }
            }

            dphidr = phidiff(currentphi, xinkey->first.get<0>()) / (_radii_all[lastgoodlayer] - _radii_all[newlayer]);
            detadr = (currenteta - xinkey->first.get<1>()) / (_radii_all[lastgoodlayer] - _radii_all[newlayer]);
            ther = (_radii_all[lastgoodlayer] - _radii_all[newlayer]) / 2;

            curvatureestimates.push_back(copysign(2 / sqrt(ther * ther + 1 / dphidr / dphidr), dphidr));

            phi_zigzag.push_back(dphidr);
            z_zigzag.push_back(detadr);

            cluskeys.push_back(xinkey->second);

            currentphi = xinkey->first.get<0>();
            currenteta = (currenteta + xinkey->first.get<1>()) / 2;

            lastgoodlayer = newlayer;
          }
        }
        if (failures > 2) continue;  //0.7 2 0.9 3

        double phi_sum = std::accumulate(phi_zigzag.begin(), phi_zigzag.end(), 0.0);
        double phi_mean = phi_sum / phi_zigzag.size();

        std::vector<double> phi_diff(phi_zigzag.size());
        std::transform(phi_zigzag.begin(), phi_zigzag.end(), phi_diff.begin(),
                       std::bind2nd(std::minus<double>(), phi_mean));
        double phi_sq_sum = std::inner_product(phi_diff.begin(), phi_diff.end(), phi_diff.begin(), 0.0);
        double phi_stdev = std::sqrt(phi_sq_sum / (phi_zigzag.size() - 1));

        double z_sum = std::accumulate(z_zigzag.begin(), z_zigzag.end(), 0.0);
        double z_mean = z_sum / z_zigzag.size();

        std::vector<double> z_diff(z_zigzag.size());
        std::transform(z_zigzag.begin(), z_zigzag.end(), z_diff.begin(),
                       std::bind2nd(std::minus<double>(), z_mean));
        double z_sq_sum = std::inner_product(z_diff.begin(), z_diff.end(), z_diff.begin(), 0.0);
        double z_stdev = std::sqrt(z_sq_sum / (z_zigzag.size() - 1));

        double curv_sum = std::accumulate(curvatureestimates.begin(), curvatureestimates.end(), 0.0);
        double curv_mean = curv_sum / curvatureestimates.size();

        std::vector<double> curv_diff(curvatureestimates.size());
        std::transform(curvatureestimates.begin(), curvatureestimates.end(), curv_diff.begin(),
                       std::bind2nd(std::minus<double>(), curv_mean));
        double curv_sq_sum = std::inner_product(curv_diff.begin(), curv_diff.end(), curv_diff.begin(), 0.0);
        double curv_stdev = std::sqrt(curv_sq_sum / (curvatureestimates.size() - 1));

        const double BQ = 0.01 * 1.4 * 0.299792458;
        double pt = BQ / abs(curv_mean);
        double pterror = BQ * curv_stdev / (curv_mean * curv_mean);
*/
        //    pt:z:dz:phi:dphi:c:dc
    // Fill NT with track parameters
    float StartPhi = atan(y0/x0);
    float StartEta = -log(tan(atan(z0/sqrt(x0*x0+y0*y0))));
    float track_pt = 1./(trackSeed.GetQPt());
    float track_pterr = trackSeed.GetErr2QPt();
    float track_z = trackSeed.GetZ();
    float track_zerr = trackSeed.GetErr2Z();
    float track_phi = 0.;
    float track_phierr = 0.001;
    float track_curvature = trackSeed.GetKappa(Bz);
    float track_curverr = 0.001;
    NT->Fill(track_pt, track_pterr, track_z, track_zerr, track_phi, track_phierr, track_curvature, track_curverr, trackKeyChain->size());
    SvtxTrack_v1 track;
    track.set_id(numberofseeds);
    for (unsigned int j = 0; j < trackKeyChain->size(); j++)
    {
      track.insert_cluster_key(trackKeyChain->at(j));
    }
    track.set_ndf(trackSeed.GetNDF());
    if(trackSeed.GetQPt()<0) track.set_charge(-1);
    else track.set_charge(1);
    TrkrCluster *cl = _cluster_map->findCluster(0);
    track.set_x(_vertex->get_x());  //track.set_x(cl->getX());
    track.set_y(_vertex->get_y());  //track.set_y(cl->getY());
    track.set_z(_vertex->get_z());  //track.set_z(cl->getZ());
    track.set_px(track_pt * cos(StartPhi));
    track.set_py(track_pt * sin(StartPhi));
    track.set_pz(track_pt / tan(2 * atan(exp(-StartEta))));
    track.set_error(0, 0, cl->getError(0, 0));
    track.set_error(0, 1, cl->getError(0, 1));
    track.set_error(0, 2, cl->getError(0, 2));
    track.set_error(1, 1, cl->getError(1, 1));
    track.set_error(1, 2, cl->getError(1, 2));
    track.set_error(2, 2, cl->getError(2, 2));
    track.set_error(3, 3, track_pterr * track_pterr * cos(StartPhi) * cos(StartPhi));
    track.set_error(4, 4, track_pterr * track_pterr * sin(StartPhi) * sin(StartPhi));
    track.set_error(5, 5, track_pterr * track_pterr / tan(2 * atan(exp(-StartEta))) / tan(2 * atan(exp(-StartEta))));
    _track_map->insert(&track);
    numberofseeds++;
  }
//t_seed->stop();
//cout << "number of seeds " << numberofseeds << endl;
//cout << "seeding time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
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
