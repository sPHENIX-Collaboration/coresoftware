/*!
 *  \file PHRTreeSeeding.C
 *  \brief Progressive pattern recgnition based on GenFit Kalman filter
 *  \detail using Alan Dion's HelixHough for seeding, GenFit Kalman filter to do track propagation
 *  \author Christof Roland & Haiwang Yu
 */

//begin

#include "PHRTreeSeeding.h"

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

// forward declarations
class PHCompositeNode;

//end

typedef bg::model::point<float, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<point, TrkrDefs::cluskey> pointKey;

#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

using namespace std;
//using namespace ROOT::Minuit2;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

vector<TrkrCluster *> clusterpoints;
SvtxVertex *_vertex;
PHRTreeSeeding::PHRTreeSeeding(
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

int PHRTreeSeeding::InitializeGeometry(PHCompositeNode *topNode)
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

bool comparing(TrkrCluster *tc1, TrkrCluster *tc2)
{
  //TrkrCluster tca = *tc1;
  //TrkrCluster tcb = *tc2;
  TrkrDefs::cluskey ck1 = tc1->getClusKey();
  TrkrDefs::cluskey ck2 = tc2->getClusKey();
  if (TrkrDefs::getLayer(ck1) != TrkrDefs::getLayer(ck2)) return TrkrDefs::getLayer(ck1) < TrkrDefs::getLayer(ck2);
  return ck1 < ck2;
  //return ck1/16<ck2/16;
}

bool largeintersection(const vector<TrkrCluster *> &v1,
                       const vector<TrkrCluster *> &v2)
{
  std::vector<TrkrCluster *> v3;
  vector<TrkrCluster *> va(v1);
  vector<TrkrCluster *> vb(v2);
  //std::sort(v1.begin(), v1.end());
  //std::sort(v2.begin(), v2.end());
  std::set_intersection(va.begin(), va.end(),
                        vb.begin(), vb.end(),
                        back_inserter(v3), comparing);
  //cout << "Got an intersection of size " << v3.size() << endl;
  if (v3.size() < 3) return false;
  //TrkrCluster *temp1 = new TrkrCluster;
  //*temp1 = *va[0];

  /*TrkrCluster *temp2 = vb[0];
  unsigned int ii = va.size();
  unsigned int jj = vb.size();
  TrkrCluster *temp3 = va[ii-1];
  TrkrCluster *temp4 = vb[jj-1];
  TrkrDefs::cluskey ck1 = (va[0])->getClusKey();
  TrkrDefs::cluskey ck2 = temp2->getClusKey();
  TrkrDefs::cluskey ck3 = temp3->getClusKey();
  TrkrDefs::cluskey ck4 = temp4->getClusKey();*/
  /*bool to_ret = (TrkrDefs::getLayer(ck1)-TrkrDefs::getLayer(ck2))*(TrkrDefs::getLayer(ck3)-TrkrDefs::getLayer(ck4))>0;
  cout << to_ret << endl;
  return to_ret;*/
  return (TrkrDefs::getLayer(va[0]->getClusKey()) - TrkrDefs::getLayer(vb[0]->getClusKey())) * (TrkrDefs::getLayer(va[va.size() - 1]->getClusKey()) - TrkrDefs::getLayer(vb[vb.size() - 1]->getClusKey())) > 0 && max(TrkrDefs::getLayer(va[va.size() - 1]->getClusKey()) - TrkrDefs::getLayer(vb[0]->getClusKey()), TrkrDefs::getLayer(vb[vb.size() - 1]->getClusKey()) - TrkrDefs::getLayer(va[0]->getClusKey())) + 1 == (int) (va.size() + vb.size() - v3.size());
}

vector<TrkrCluster *> yunion(const vector<TrkrCluster *> &v1, const vector<TrkrCluster *> &v2)
{
  std::vector<TrkrCluster *> v3;

  //std::sort(v1.begin(), v1.end());
  //std::sort(v2.begin(), v2.end());

  std::set_union(v1.begin(), v1.end(),
                 v2.begin(), v2.end(),
                 back_inserter(v3), comparing);
  return v3;
}

bool issuperset(const vector<TrkrCluster *> &larger,
                const vector<TrkrCluster *> &smaller)
{
  for (unsigned int i = 0; i < smaller.size(); i++)
  {
    if (!binary_search(larger.begin(), larger.end(), smaller.at(i), comparing)) return false;
  }
  return true;
}

void changetounion(vector<TrkrCluster *> &v1,
                   const vector<TrkrCluster *> &v2)
{
  vector<TrkrCluster *> temp(yunion(v1, v2));
  v1 = temp;
}

double RosenBrock(const double *xx)
{
  /* const Double_t x = xx[0];
  const Double_t y = xx[1];
  const Double_t tmp1 = y-x*x;
  const Double_t tmp2 = 1-x;
  return 100*tmp1*tmp1+tmp2*tmp2;*/
  //printf("arguments: %f %f %f %f %f", xx[0],xx[1],xx[2],xx[3],xx[4]);
  double chi2 = 0.;
  for (unsigned int i = 0; i < clusterpoints.size(); i++)
  {
    TrkrCluster *tc = clusterpoints.at(i);
    double qz = xx[i + 5] - clusterpoints.at(0)->getZ();
    double predictedx = xx[0] - sin(xx[2]) / xx[4] + cos(xx[4] / tan(xx[3]) * qz + xx[2] - M_PI / 2) / xx[4];
    double predictedy = xx[1] + cos(xx[2]) / xx[4] + sin(xx[4] / tan(xx[3]) * qz + xx[2] - M_PI / 2) / xx[4];
    double xresid = tc->getX() - predictedx;
    double yresid = tc->getY() - predictedy;
    double zresid = tc->getZ() - xx[i + 5];
    //cout << "zresid=" << zresid << endl;
    double contribution = xresid * xresid / tc->getError(0, 0) + yresid * yresid / tc->getError(1, 1) + zresid * zresid / tc->getError(2, 2);
    //cout << "contribution of " << contribution << endl;
    chi2 += contribution;
  }
  double qz = xx[clusterpoints.size() + 5] - _vertex->get_z();
  double predictedx = xx[0] - sin(xx[2]) / xx[4] + cos(xx[4] / tan(xx[3]) * qz + xx[2] - M_PI / 2) / xx[4];
  double predictedy = xx[1] + cos(xx[2]) / xx[4] + sin(xx[4] / tan(xx[3]) * qz + xx[2] - M_PI / 2) / xx[4];
  double xresid = _vertex->get_x() - predictedx;
  double yresid = _vertex->get_y() - predictedy;
  double zresid = _vertex->get_z() - xx[clusterpoints.size() + 5];
  //cout << "zresid=" << zresid << endl;
  double contribution = xresid * xresid / _vertex->get_error(0, 0) + yresid * yresid / _vertex->get_error(1, 1) + zresid * zresid / _vertex->get_error(2, 2);
  //cout << "contribution of " << contribution << endl;
  chi2 += contribution;
  return chi2;
}

int PHRTreeSeeding::GetNodes(PHCompositeNode *topNode)
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

double PHRTreeSeeding::phiadd(double phi1, double phi2)
{
  double s = phi1 + phi2;
  if (s > 2 * M_PI)
    return s - 2 * M_PI;
  else if (s < 0)
    return s + 2 * M_PI;
  else
    return s;
}

double PHRTreeSeeding::phidiff(double phi1, double phi2)
{
  double d = phi1 - phi2;
  if (d > M_PI)
    return d - 2 * M_PI;
  else if (d < -M_PI)
    return d + 2 * M_PI;
  else
    return d;
}

void PHRTreeSeeding::QueryTree(const bgi::rtree<pointKey, bgi::quadratic<16>> &rtree, double phimin, double etamin, double lmin, double phimax, double etamax, double lmax, std::vector<pointKey> &returned_values)
{
  rtree.query(bgi::intersects(box(point(phimin, etamin, lmin), point(phimax, etamax, lmax))), std::back_inserter(returned_values));
  if (phimin < 0) rtree.query(bgi::intersects(box(point(2 * M_PI + phimin, etamin, lmin), point(2 * M_PI, etamax, lmax))), std::back_inserter(returned_values));
  if (phimax > 2 * M_PI) rtree.query(bgi::intersects(box(point(0, etamin, lmin), point(phimax - 2 * M_PI, etamax, lmax))), std::back_inserter(returned_values));
}

/*double PHRTreeSeeding::pointKeyToTuple(pointKey *pK)
{
  TVector3 *v;
  v->SetPtEtaPhi(_radii_all[(int)(pK->first.get<2>+0.4)], pK->first.get<1>, pK->first.get<0>);
  return make_tuple(v.X(), v.Y(), v.Z());
  }*/

double /*PHRTreeSeeding::*/ costfunction(const double *xx)
// xcenter, ycenter, radius, frequency, phase
{
  /*double cost = 0;
  //should add a term for the vertex
  for(vector<TrkrCluster*>::iterator i = clusterpoints.begin(); i != clusterpoints.end(); i++)
  {
    TrkrCluster* tc = *i;
    double qa = tc->getError(0,0);
    double qb = tc->getError(0,1);
    double qc = tc->getError(1,0);
    double qd = tc->getError(1,1);
    double qz = tc->getZ();
    double qx = tc->getX()-xx[0]-xx[2]*cos(xx[3]*qz+xx[4]);
    double qy = tc->getY()-xx[1]-xx[2]*sin(xx[3]*qz+xx[4]);

    cost += (qd*qx*qx-(qb+qc)*qx*qy+qa*qy*qy)/(qa*qd-qb*qc);
      //cost += pow(get<0>(clusterpoints.at(i))-xx[0]-xx[2]*cos(xx[3]*get<2>(clusterpoints.at(i))+xx[4]),2)+pow(get<1>(clusterpoints.at(i))-xx[1]-xx[2]*sin(xx[3]*get<2>(clusterpoints.at(i))+xx[4]),2)
  }
  return cost;*/
      return xx[0] + xx[1] + xx[2] + xx[3] + xx[4];
}

double PHRTreeSeeding::chisq(const double *xx)
{
  double chi2 = 0.;
  for (vector<TrkrCluster *>::iterator i = clusterpoints.begin(); i != clusterpoints.end(); i++)
  {
    TrkrCluster *tc = *i;
    double qa = tc->getError(0, 0);
    double qb = tc->getError(0, 1);
    double qc = tc->getError(1, 0);
    double qd = tc->getError(1, 1);
    double qz = tc->getZ() - clusterpoints.at(0)->getZ();
    double predictedx = xx[0] - sin(xx[2]) / xx[4] + cos(xx[4] / tan(xx[3]) * qz + xx[2] - M_PI / 2) / xx[4];
    double predictedy = xx[1] + cos(xx[2]) / xx[4] + sin(xx[4] / tan(xx[3]) * qz + xx[2] - M_PI / 2) / xx[4];
    double qx = tc->getX() - predictedx;
    double qy = tc->getY() - predictedy;

    chi2 += (qd * qx * qx - (qb + qc) * qx * qy + qa * qy * qy) / (qa * qd - qb * qc);
  }

  return chi2;
}

void PHRTreeSeeding::FillTree()
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

int PHRTreeSeeding::Process(PHCompositeNode *topNode)
{
  TFile fpara("rtree_para.root", "RECREATE");
  TNtuple *NT = new TNtuple("NT", "NT", "pt:dpt:z:dz:phi:dphi:c:dc:nhit");
  _vertex = _vertex_map->get(0);

  //for different purpose
  phisr = 0.005 * _phi_scale;  // *2
  etasr = 0.0035 * _z_scale;   // *2;
  /* 0.7 version
     phist = 0.001*1;
     etast = 0.003*1;
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
  cout << " entries in tree: " << _rtree.size() << endl;

  for (unsigned int iteration = 0; iteration < 1; iteration++)
  {
    if (iteration == 1) _start_layer -= 7;
    vector<pointKey> StartLayerClusters;
    _rtree.query(bgi::intersects(
                     box(point(0, -3, ((double) _start_layer - 0.5)),
                         point(2 * M_PI, 3, ((double) _start_layer + 0.5)))),
                 std::back_inserter(StartLayerClusters));

    for (vector<pointKey>::iterator StartCluster = StartLayerClusters.begin(); StartCluster != StartLayerClusters.end(); StartCluster++)
    {
      double StartPhi = StartCluster->first.get<0>();
      double StartEta = StartCluster->first.get<1>();

      vector<pointKey> SecondLayerClusters;
      QueryTree(_rtree,
                StartPhi - phisr,
                StartEta - etasr,
                (double) _start_layer - 1.5,
                StartPhi + phisr,
                StartEta + etasr,
                (double) _start_layer - 0.5,
                SecondLayerClusters);
      cout << " entries in second layer: " << SecondLayerClusters.size() << endl;

      for (vector<pointKey>::iterator SecondCluster = SecondLayerClusters.begin(); SecondCluster != SecondLayerClusters.end(); SecondCluster++)
      {
        double currentphi = SecondCluster->first.get<0>();
        double currenteta = SecondCluster->first.get<1>();
        int lastgoodlayer = _start_layer - 1;

        int failures = 0;

        double dphidr = phidiff(StartPhi, currentphi) / (_radii_all[_start_layer] - _radii_all[_start_layer - 1]);
        double detadr = (StartEta - currenteta) / (_radii_all[_start_layer] - _radii_all[_start_layer - 1]);

        double ther = (_radii_all[_start_layer] + _radii_all[_start_layer - 1]) / 2;

        vector<double> curvatureestimates;
        curvatureestimates.push_back(copysign(2 / sqrt(ther * ther + 1 / dphidr / dphidr), dphidr));
        vector<double> phi_zigzag;
        phi_zigzag.push_back(dphidr);
        vector<double> z_zigzag;
        z_zigzag.push_back(detadr);

        vector<TrkrDefs::cluskey> cluskeys;
        cluskeys.push_back(StartCluster->second);
        cluskeys.push_back(SecondCluster->second);
        cout << " phi 1: " << StartPhi
             << " phi 2: " << currentphi
             << " dphi: " << dphidr
             << " eta 1: " << StartEta
             << " eta 2: " << currenteta
             << " deta: " << detadr
             << endl;

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
            for (std::vector<pointKey>::iterator it = newlayer_clusters.begin(); it != newlayer_clusters.end(); it++)
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

        //    pt:z:dz:phi:dphi:c:dc
        NT->Fill(pt, pterror, z_mean, z_stdev, phi_mean, phi_stdev, curv_mean, curv_stdev, cluskeys.size());
        SvtxTrack_v1 track;
        track.set_id(numberofseeds);

        for (unsigned int j = 0; j < cluskeys.size(); j++)
        {
          track.insert_cluster_key(cluskeys.at(j));
        }
        track.set_ndf(2 * cluskeys.size() - 5);
        short int helicity = 1;
        if (StartPhi * curv_mean < 0) helicity -= 2;
        track.set_charge(-helicity);

        TrkrCluster *cl = _cluster_map->findCluster(StartCluster->second);
        track.set_x(_vertex->get_x());  //track.set_x(cl->getX());
        track.set_y(_vertex->get_y());  //track.set_y(cl->getY());
        track.set_z(_vertex->get_z());  //track.set_z(cl->getZ());
        track.set_px(pt * cos(StartPhi));
        track.set_py(pt * sin(StartPhi));
        track.set_pz(pt / tan(2 * atan(exp(-StartEta))));
        track.set_error(0, 0, cl->getError(0, 0));
        track.set_error(0, 1, cl->getError(0, 1));
        track.set_error(0, 2, cl->getError(0, 2));
        track.set_error(1, 1, cl->getError(1, 1));
        track.set_error(1, 2, cl->getError(1, 2));
        track.set_error(2, 2, cl->getError(2, 2));
        track.set_error(3, 3, pterror * pterror * cos(StartPhi) * cos(StartPhi));
        track.set_error(4, 4, pterror * pterror * sin(StartPhi) * sin(StartPhi));
        track.set_error(5, 5, pterror * pterror / tan(2 * atan(exp(-StartEta))) / tan(2 * atan(exp(-StartEta))));
        _track_map->insert(&track);

        numberofseeds++;
      }
    }
  }
  t_seed->stop();
  cout << "number of seeds " << numberofseeds << endl;
  cout << "seeding time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;

  /*

  */

  fpara.cd();
  NT->Write();
  fpara.Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHRTreeSeeding::Setup(PHCompositeNode *topNode)
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

int PHRTreeSeeding::End()
{
  cout << "Called End " << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
