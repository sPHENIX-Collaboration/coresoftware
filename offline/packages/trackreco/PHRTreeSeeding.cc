
 /*!
 *  \file PHRTreeSeeding.C
 *  \brief Progressive pattern recgnition based on GenFit Kalman filter
 *  \detail using Alan Dion's HelixHough for seeding, GenFit Kalman filter to do track propagation
 *  \author Christof Roland & Haiwang Yu
 */



//begin

//#ifndef G4HOUGH_PHG4KALMANPATREC_H
//#define G4HOUGH_PHG4KALMANPATREC_H

#include <fun4all/SubsysReco.h>

// Helix Hough + Eigen includes (hidden from rootcint)
#if !defined(__CINT__) || defined(__CLING__)
#include <HelixHough/SimpleHit3D.h>
#include <HelixHough/SimpleTrack3D.h>
#include <HelixHough/VertexFinder.h>
#include <Eigen/Core>                              // for Matrix
#endif



#if !defined(__CINT__) || defined(__CLING__)
//BOOST for combi seeding
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>

#include <boost/geometry/index/rtree.hpp>
#endif

// standard includes
#include <cfloat>
#include <iostream>                                // for operator<<, basic_...
#include <list>
#include <map>
#include <memory>
#include <set>                                     // for set
#include <string>                                  // for string
#include <utility>                                 // for pair
#include <vector>

// forward declarations
class BbcVertexMap;

class PHCompositeNode;

class PHG4CellContainer;
class PHG4CylinderGeomContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;
class PHTimer;

class sPHENIXSeedFinder;

class SvtxClusterMap;
class SvtxCluster;
class SvtxTrackMap;
class SvtxTrack;
class SvtxTrackState;
class SvtxVertexMap;
class SvtxHitMap;

class TNtuple;
class TFile;

namespace PHGenFit
{
class Fitter;
class Track;
class Measurement;
} /* namespace PHGenFit */
//end



#include "PHRTreeSeeding.h"

#include "AssocInfoContainer.h"                         // for AssocInfoCont...

// Helix Hough includes
#include <HelixHough/HelixKalmanState.h>                // for HelixKalmanState
#include <HelixHough/HelixRange.h>
#include <HelixHough/SimpleHit3D.h>
#include <HelixHough/SimpleTrack3D.h>
#include <HelixHough/sPHENIXSeedFinder.h>               // for sPHENIXSeedFi...
#include <HelixHough/VertexFinder.h>


// trackbase_historic includes
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex_v1.h>

#include <trackbase/TrkrCluster.h>                      // for TrkrCluster
#include <trackbase/TrkrDefs.h>                         // for getLayer, clu...
#include <trackbase/TrkrClusterContainer.h>

// sPHENIX Geant4 includes
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

// sPHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHTimer.h>                              // for PHTimer
#include <phool/getClass.h>

#include <Eigen/Core>                  // for Matrix
#include <Eigen/Dense>

//Helix FCN includes
//#include "Minuit2-5.34.14/test/MnSim/GaussFcn.h"
//#include "Minuit2-5.34.14/test/MnSim/HelixFcn.h"
//#include "Minuit2-5.34.14/test/MnSim/HelixFcn.cpp"
//#include "HelixFcn.cpp"






//#include "Minuit2-5.34.14/inc/Minuit2/FunctionMinimum.h"
//#include "Minuit2-5.34.14/inc/Minuit2/MnUserParameters.h"
//#include "Minuit2-5.34.14/inc/Minuit2/MnUserParameterState.h"
//#include "Minuit2-5.34.14/inc/Minuit2/MinimumPrint.h"
//#include "Minuit2-5.34.14/inc/Minuit2/MnMigrad.h"
//#include "Minuit2-5.34.14/inc/Minuit2/MnMinos.h"
//#include "Minuit2-5.34.14/inc/Minuit2/MnContours.h"
//#include "Minuit2-5.34.14/inc/Minuit2/MnPlot.h"
//#include "Minuit2-5.34.14/inc/Minuit2/MnStrategy.h"
// #include "Minuit2/FunctionMinimum.h"
// #include "Minuit2/MnUserParameterState.h"
// #include "Minuit2/MnPrint.h"
// #include "Minuit2/MnMigrad.h"
// #include "Minuit2/MnMinos.h"
// #include "Minuit2/MnContours.h"
// #include "Minuit2/MnPlot.h"
// #include "Minuit2/MinosError.h"
// #include "Minuit2/ContoursError.h"



//ROOT includes for debugging
#include <TFile.h>
#include <TNtuple.h>
#include <TAxis.h>
#include <TGraph.h>
//BOOST for combi seeding
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

typedef bg::model::point<float, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<point, TrkrDefs::cluskey> pointKey;



// standard includes
#include <TH1.h>
#include <stdio.h>      /* printf */
#include <math.h>       /* copysign */
#include <algorithm>
#include <climits>                                     // for UINT_MAX
#include <cmath>
#include <iostream>
#include <memory>
#include <set>                                          // for set
#include <tuple>
#include <utility>   
                                   // for pair, make_pair
//#include <assert.h>

//#include "Math/Minimizer.h"
//#include "Math/Factory.h"
//#include "Math/Functor.h"
//#include "Minuit2/Minuit2Minimizer.h"
//#include "Math/Functor.h"




#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp


using namespace std;
//using namespace ROOT::Minuit2;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

//typedef uint64_t cluskey;

vector<TrkrCluster*> clusterpoints;
SvtxVertex *_vertex;
PHRTreeSeeding::PHRTreeSeeding(
			       const string& name,
			       unsigned int nlayers_maps,
			       unsigned int nlayers_intt,
			       unsigned int nlayers_tpc,
			       unsigned int nlayers_seeding,
			       unsigned int min_nlayers_seeding)//:
//_g4clusters(nullptr)
{
  
}



/*
int PHRTreeSeeding::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------
  cout << "called GetNodes" << endl;

  _g4clusters = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
  if (!_g4clusters)
  {
    cerr << PHWHERE << " ERROR: Can't find node SvtxClusterMap" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return  Fun4AllReturnCodes::EVENT_OK;
}
*/

int PHRTreeSeeding::InitializeGeometry(PHCompositeNode* topNode)
{


  PHG4CylinderCellGeomContainer* cellgeos = findNode::getClass<
      PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  PHG4CylinderGeomContainer* laddergeos = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  PHG4CylinderGeomContainer* mapsladdergeos = findNode::getClass<
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
      PHG4CylinderCellGeom* geo = miter->second;
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
      PHG4CylinderGeom* geo = miter->second;
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
      PHG4CylinderGeom* geo = miter->second;

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

bool comparing(TrkrCluster* tc1, TrkrCluster* tc2)
{
  //TrkrCluster tca = *tc1;
  //TrkrCluster tcb = *tc2;
  TrkrDefs::cluskey ck1 = tc1->getClusKey();
  TrkrDefs::cluskey ck2 = tc2->getClusKey();
  if(TrkrDefs::getLayer(ck1)!=TrkrDefs::getLayer(ck2)) return TrkrDefs::getLayer(ck1)<TrkrDefs::getLayer(ck2);
  return ck1<ck2;
  //return ck1/16<ck2/16;
}


bool largeintersection(const vector<TrkrCluster*> &v1,
		       const vector<TrkrCluster*> &v2){
  std::vector<TrkrCluster*> v3;
  vector<TrkrCluster*> va(v1);
  vector<TrkrCluster*> vb(v2);
  //std::sort(v1.begin(), v1.end());
  //std::sort(v2.begin(), v2.end());
  std::set_intersection(va.begin(),va.end(),
			vb.begin(),vb.end(),
			back_inserter(v3),comparing);
  //cout << "Got an intersection of size " << v3.size() << endl;
  if(v3.size()<3) return false;
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
    return (TrkrDefs::getLayer(va[0]->getClusKey())-TrkrDefs::getLayer(vb[0]->getClusKey()))*(TrkrDefs::getLayer(va[va.size()-1]->getClusKey())-TrkrDefs::getLayer(vb[vb.size()-1]->getClusKey()))>0 && max(TrkrDefs::getLayer(va[va.size()-1]->getClusKey())-TrkrDefs::getLayer(vb[0]->getClusKey()),TrkrDefs::getLayer(vb[vb.size()-1]->getClusKey())-TrkrDefs::getLayer(va[0]->getClusKey()))+1==(int)(va.size()+vb.size()-v3.size());
}

vector<TrkrCluster*> yunion(const vector<TrkrCluster*> &v1, const vector<TrkrCluster*> &v2){
  std::vector<TrkrCluster*> v3;
  
  //std::sort(v1.begin(), v1.end());
  //std::sort(v2.begin(), v2.end());

  std::set_union(v1.begin(),v1.end(),
		 v2.begin(),v2.end(),
		 back_inserter(v3),comparing);
  return v3;
}

bool issuperset(const vector<TrkrCluster*> &larger,
		const vector<TrkrCluster*> &smaller){
  for(unsigned int i=0; i<smaller.size();i++){
    if(!binary_search(larger.begin(),larger.end(),smaller.at(i),comparing)) return false;
  }
  return true;
}



void changetounion(vector<TrkrCluster*> &v1,
		   const vector<TrkrCluster*> &v2){
  vector<TrkrCluster*> temp(yunion(v1,v2));
  v1 = temp;
}






double RosenBrock(const double *xx )
{
  /* const Double_t x = xx[0];
  const Double_t y = xx[1];
  const Double_t tmp1 = y-x*x;
  const Double_t tmp2 = 1-x;
  return 100*tmp1*tmp1+tmp2*tmp2;*/
  //printf("arguments: %f %f %f %f %f", xx[0],xx[1],xx[2],xx[3],xx[4]);
  double chi2 = 0.;
  for(unsigned int i=0;i<clusterpoints.size();i++){
    TrkrCluster* tc = clusterpoints.at(i);
    double qz = xx[i+5]-clusterpoints.at(0)->getZ();
    double predictedx = xx[0]-sin(xx[2])/xx[4]+cos(xx[4]/tan(xx[3])*qz+xx[2]-M_PI/2)/xx[4];
    double predictedy = xx[1]+cos(xx[2])/xx[4]+sin(xx[4]/tan(xx[3])*qz+xx[2]-M_PI/2)/xx[4];
    double xresid = tc->getX()-predictedx;
    double yresid = tc->getY()-predictedy;
    double zresid = tc->getZ()-xx[i+5];
    //cout << "zresid=" << zresid << endl;
    double contribution = xresid*xresid/tc->getError(0,0)+yresid*yresid/tc->getError(1,1)+zresid*zresid/tc->getError(2,2);
    //cout << "contribution of " << contribution << endl;
    chi2 += contribution;
  }
  double qz = xx[clusterpoints.size()+5]-_vertex->get_z();
  double predictedx = xx[0]-sin(xx[2])/xx[4]+cos(xx[4]/tan(xx[3])*qz+xx[2]-M_PI/2)/xx[4];
  double predictedy = xx[1]+cos(xx[2])/xx[4]+sin(xx[4]/tan(xx[3])*qz+xx[2]-M_PI/2)/xx[4];
  double xresid = _vertex->get_x()-predictedx;
  double yresid = _vertex->get_y()-predictedy;
  double zresid = _vertex->get_z()-xx[clusterpoints.size()+5];
    //cout << "zresid=" << zresid << endl;
  double contribution = xresid*xresid/_vertex->get_error(0,0)+yresid*yresid/_vertex->get_error(1,1)+zresid*zresid/_vertex->get_error(2,2);
    //cout << "contribution of " << contribution << endl;
    chi2 += contribution;
  return chi2;
}

int PHRTreeSeeding::GetNodes(PHCompositeNode* topNode)
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

double PHRTreeSeeding::phiadd(double phi1, double phi2){
  double s = phi1+phi2;
  if(s>2*M_PI) return s-2*M_PI;
  else if(s<0) return s+2*M_PI;
  else return s;
}


double PHRTreeSeeding::phidiff(double phi1, double phi2){
  double d = phi1-phi2;
  if(d>M_PI) return d-2*M_PI;
  else if(d<-M_PI) return d+2*M_PI;
  else return d;
}


void PHRTreeSeeding::wquery(const bgi::rtree<pointKey, bgi::quadratic<16>> &rtree, double phimin, double etamin, double lmin, double phimax, double etamax, double lmax,std::vector<pointKey> &returned_values){
  rtree.query(bgi::intersects(box(point(phimin,etamin,lmin),point(phimax,etamax,lmax))),std::back_inserter(returned_values));
  if(phimin<0) rtree.query(bgi::intersects(box(point(2*M_PI+phimin,etamin,lmin),point(2*M_PI,etamax,lmax))),std::back_inserter(returned_values));
  if(phimax>2*M_PI) rtree.query(bgi::intersects(box(point(0,etamin,lmin),point(phimax-2*M_PI,etamax,lmax))),std::back_inserter(returned_values));
}


/*double PHRTreeSeeding::pointKeyToTuple(pointKey *pK)
{
  TVector3 *v;
  v->SetPtEtaPhi(_radii_all[(int)(pK->first.get<2>+0.4)], pK->first.get<1>, pK->first.get<0>);
  return make_tuple(v.X(), v.Y(), v.Z());
  }*/

double /*PHRTreeSeeding::*/costfunction(const double *xx)
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
  return cost;*/return xx[0]+xx[1]+xx[2]+xx[3]+xx[4];
}

double PHRTreeSeeding::chisq(const double *xx)
{
  double chi2 = 0.;
  for(vector<TrkrCluster*>::iterator i = clusterpoints.begin(); i != clusterpoints.end(); i++) {
    TrkrCluster* tc = *i;
    double qa = tc->getError(0,0);
    double qb = tc->getError(0,1);
    double qc = tc->getError(1,0);
    double qd = tc->getError(1,1);
    double qz = tc->getZ()-clusterpoints.at(0)->getZ();
    double predictedx = xx[0]-sin(xx[2])/xx[4]+cos(xx[4]/tan(xx[3])*qz+xx[2]-M_PI/2)/xx[4];
    double predictedy = xx[1]+cos(xx[2])/xx[4]+sin(xx[4]/tan(xx[3])*qz+xx[2]-M_PI/2)/xx[4];
    double qx = tc->getX()-predictedx;
    double qy = tc->getY()-predictedy;

    chi2 += (qd*qx*qx-(qb+qc)*qx*qy+qa*qy*qy)/(qa*qd-qb*qc);
  }

  return chi2;
}

/*int PHRTreeSeeding::Process()
{
  return  Fun4AllReturnCodes::EVENT_OK;
}
*/

int PHRTreeSeeding::Process(PHCompositeNode *topNode)
{

  _vertex = _vertex_map->get(0);



  /*phisr = 0.1;
  etasr = 0.015;
  phist = 0.04;
  etast = 0.01;
  phixt = 0.03;
  etaxt = 0.008;*/
  //phisr = 0.015;
  //etasr = 0.005;
  //phist = 0.01;
  //etast = 0.008;
  //phixt = 0.008;
  //etaxt = 0.005;


  //for different purpose
  phisr = 0.005*2;
  etasr = 0.0035*2;
  phist = 0.001*1;
  etast = 0.003*1;

  bgi::rtree<pointKey, bgi::quadratic<16> > rtree;
  PHTimer* t_fill = new PHTimer("t_fill");
  t_fill->stop();
  //  int num = 0;
  int nlayer[60];
  for (int j = 0; j < 60; j++) nlayer[j] = 0;
  TrkrClusterContainer::ConstRange clusrange = _cluster_map->getClusters();
  TGraph *g = new TGraph();
  g->GetXaxis()->SetTitle("r");
  g->GetYaxis()->SetTitle("_radii_all[layer]");
  //int laod = 0;
  //TH1D *phierrors = new TH1D("phierrors","phierrors",100,0,0.0016);
  //TH1D *etaerrors = new TH1D("etaerrors","etaerrors",100,0,0.004);
  for(TrkrClusterContainer::ConstIterator iter = clusrange.first; iter != clusrange.second; ++iter)  
  {
    /*if (_hit_used_map[iter->first] != 0)
    {
      continue;
      }*/
    TrkrCluster* cluster = iter->second;
    TrkrDefs::cluskey ckey = iter->first;
    unsigned int layer = TrkrDefs::getLayer(ckey);
    if(layer<44) continue;
    /*if (layer < _nlayers_maps + _nlayers_intt)
      continue;*/
    TVector3 vec(cluster->getPosition(0)-_vertex->get_x(), cluster->getPosition(1)-_vertex->get_y(), cluster->getPosition(2)-_vertex->get_z());


    /*double CVxx = cluster->getError(0,0);
    double CVxy = cluster->getError(0,1);
    double CVxz = cluster->getError(0,2);
    double CVyy = cluster->getError(1,1);
    double CVyz = cluster->getError(1,2);
    double CVzz = cluster->getError(2,2);*/
    
    //phierrors->Fill(sqrt((CVyy*vec.X()*vec.X()-2*CVxy*vec.X()*vec.Y()+CVxx*vec.Y()*vec.Y())/((vec.X()*vec.X()+vec.Y()*vec.Y())*(vec.X()*vec.X()+vec.Y()*vec.Y()))));
    //etaerrors->Fill(sqrt((CVzz*(vec.X()*vec.X()+vec.Y()*vec.Y())*(vec.X()*vec.X()+vec.Y()*vec.Y())+vec.Z()*(-2*(CVxz*vec.X()+CVyz*vec.Y())*(vec.X()*vec.X()+vec.Y()*vec.Y())+CVxx*vec.X()*vec.X()*vec.Z()+CVyy*vec.Y()*vec.Y()*vec.Z()+2*CVxy*vec.X()*vec.Y()*vec.Z()))/((vec.X()*vec.X()+vec.Y()*vec.Y())*(vec.X()*vec.X()+vec.Y()*vec.Y())*(vec.X()*vec.X()+vec.Y()*vec.Y()+vec.Z()*vec.Z()))));


    //floats or doubles?
    //cout << vec.X() << "," << vec.Y() << "," << vec.Z() << endl;
    double clus_phi = vec.Phi();
    clus_phi -= 2*M_PI*floor(clus_phi/(2*M_PI));
    double clus_eta = vec.Eta();
    double  clus_l =layer;// _radii_all[layer];
    vector<pointKey> testduplicate;
    wquery(rtree,clus_phi-0.00001,clus_eta-0.00001,layer-0.5,clus_phi+0.00001,clus_eta+0.00001,layer+0.5,testduplicate);
    if(!testduplicate.empty()) continue;
    nlayer[layer]++;
    t_fill->restart();
    rtree.insert(std::make_pair(point(clus_phi, clus_eta, clus_l), ckey));
    t_fill->stop();
    //double myx = sqrt(cluster->getPosition(0)*cluster->getPosition(0)+cluster->getPosition(1)*cluster->getPosition(1)+0*cluster->getPosition(2)*cluster->getPosition(2));
    //double myy = _radii_all[layer];
    //g->SetPoint(laod,myx,myy-myx);
    //laod++;
  }
  //g->Draw("*");
  TFile *f = new TFile("/sphenix/user/rkao/RBAREunModular/test/dontusephietaerror.root", "RECREATE");
  f->cd();
  //g->Write("heeeh");
  //phierrors->Write("phierrors");
  //etaerrors->Write("etaerrors");


  std::cout << "fill time: " << t_fill->get_accumulated_time() / 1000. << " sec" << std::endl;
  
  PHTimer* t_seed = new PHTimer("t_seed");
  t_seed->stop();
  t_seed->restart();


  int numberofseeds = 0;

  vector<pointKey> layer53clusters;
  rtree.query(bgi::intersects(box(point(0,-3,52.5),point(2*M_PI,3,53.5))),std::back_inserter(layer53clusters));
  for(vector<pointKey>::iterator layer53cluster = layer53clusters.begin();layer53cluster!=layer53clusters.end();layer53cluster++){
    double layer53phi = layer53cluster->first.get<0>();
    double layer53eta = layer53cluster->first.get<1>();
    //cout << "numberofseeds is now " << numberofseeds << endl;
    vector<pointKey> layer52clusters;
    wquery(rtree,layer53phi-phisr,layer53eta-etasr,51.5,layer53phi+phisr,layer53eta+etasr,52.5,layer52clusters);
    for(vector<pointKey>::iterator layer52cluster = layer52clusters.begin();layer52cluster!=layer52clusters.end();layer52cluster++){
      double currentphi = layer52cluster->first.get<0>();
      double currenteta = layer52cluster->first.get<1>();
      int lastgoodlayer = 52;
      //int newlayer = 51;
      int failures = 0;
      double dphidr = phidiff(layer53phi,currentphi)/(_radii_all[53]-_radii_all[52]);
      double ther = (_radii_all[53]+_radii_all[52])/2;
      vector<double> curvatureestimates;
      curvatureestimates.push_back(copysign(2/sqrt(ther*ther+1/dphidr/dphidr),dphidr));
      vector<TrkrDefs::cluskey> cluskeys;
      cluskeys.push_back(layer53cluster->second);
      cluskeys.push_back(layer52cluster->second);
      //while(newlayer>=44){
      for(int newlayer = 51; newlayer>=44;newlayer--){
	//cout << "newlayer is now " << newlayer << endl;
	vector<pointKey> newlayer_clusters;
	wquery(rtree,currentphi-dphidr*(_radii_all[lastgoodlayer]-_radii_all[newlayer])-phist,currenteta-etast,newlayer-0.5,currentphi-dphidr*(_radii_all[lastgoodlayer]-_radii_all[newlayer])+phist,currenteta+etast,newlayer+0.5,newlayer_clusters);
	if(newlayer_clusters.empty()){
	  failures+=1;
	  if(failures>3) break;
	}
	else{
	  double xinrecord = 100.0;
	  pointKey *xinkey = &*newlayer_clusters.begin();
	  for(std::vector<pointKey>::iterator it=newlayer_clusters.begin(); it!=newlayer_clusters.end();it++){
	    double dist = abs(phidiff(it->first.get<0>(),currentphi-dphidr*(_radii_all[lastgoodlayer]-_radii_all[newlayer])))+abs(it->first.get<1>()-currenteta);
	    if(dist < xinrecord){
	      *xinkey = *it;
	      xinrecord = dist;
	    }
	  }

	  dphidr = phidiff(currentphi,xinkey->first.get<0>())/(_radii_all[lastgoodlayer]-_radii_all[newlayer]);
	  ther = (_radii_all[lastgoodlayer]-_radii_all[newlayer])/2;
	  curvatureestimates.push_back(copysign(2/sqrt(ther*ther+1/dphidr/dphidr),dphidr));
	  cluskeys.push_back(xinkey->second);
	  currentphi = xinkey->first.get<0>();
	  currenteta = (currenteta+xinkey->first.get<1>())/2;
	  lastgoodlayer = newlayer;
	}
	//newlayer--;
      }
      if(failures>3) continue;

      double sum = std::accumulate(curvatureestimates.begin(), curvatureestimates.end(), 0.0);
      double mean = sum / curvatureestimates.size();

      std::vector<double> diff(curvatureestimates.size());
      std::transform(curvatureestimates.begin(), curvatureestimates.end(), diff.begin(),
		     std::bind2nd(std::minus<double>(), mean));
      double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
      double stdev = std::sqrt(sq_sum / (curvatureestimates.size()-1));
      
      const double BQ = 0.01*1.4*0.299792458;
      double pt = BQ/abs(mean);
      double pterror = BQ*stdev/(mean*mean);


      SvtxTrack_v1 track;
      track.set_id(numberofseeds);
      //for(unsigned int j = 0; j<clusterpoints.size();j++) track.insert_cluster_key(clusterpoints.at(j)->getClusKey());
      //cout << "cluskeys size is " << cluskeys.size() << endl;
      for(unsigned int j=0; j<cluskeys.size(); j++) {track.insert_cluster_key(cluskeys.at(j));/*printf("A cluskey: %lX", cluskeys.at(j));*/}
      //track.set_chisq(500);
      track.set_ndf(2*cluskeys.size()-5);
      short int helicity = 1;
      if(layer53eta*mean<0) helicity-=2;
      track.set_charge(-helicity);

      TrkrCluster *cl = _cluster_map->findCluster(layer53cluster->second);
      track.set_x(cl->getX());
      track.set_y(cl->getY());
      track.set_y(cl->getZ());
      track.set_px(pt*cos(layer53phi));
      track.set_py(pt*sin(layer53phi));
      track.set_pz(pt/tan(2*atan(exp(-layer53eta))));
      track.set_error(0,0,cl->getError(0,0));
      track.set_error(0,1,cl->getError(0,1));
      track.set_error(0,2,cl->getError(0,2));
      track.set_error(1,1,cl->getError(1,1));
      track.set_error(1,2,cl->getError(1,2));
      track.set_error(2,2,cl->getError(2,2));
      track.set_error(3,3,pterror*pterror*cos(layer53phi)*cos(layer53phi));
      track.set_error(4,4,pterror*pterror*sin(layer53phi)*sin(layer53phi));
      track.set_error(5,5,pterror*pterror/tan(2*atan(exp(-layer53eta)))/tan(2*atan(exp(-layer53eta))));
      //printf("Found a track with px=%f+/-%f, py=%f+/-%f, pz=%f+/-%f\n", pt*cos(layer53phi),pterror*cos(layer53phi),pt*sin(layer53phi),pterror*sin(layer53phi),pt/tan(2*atan(exp(-layer53eta))),pterror/tan(2*atan(exp(-layer53eta))));
      _track_map->insert(&track);


      numberofseeds++;
    }
  }
  t_seed->stop();
  //cout << "number of seeds " << numberofseeds << endl;
  //cout << "seeding time: " << t_seed->get_accumulated_time()/1000 << " s" << endl;

  /*for(unsigned int i=0;i<_track_map->size();i++){
    printf("Found a %d cluster track  with px=%f, py=%f, pz=%f\n", (int)_track_map->get(i)->size_cluster_keys(),_track_map->get(i)->get_px(),_track_map->get(i)->get_py(),_track_map->get(i)->get_pz());
    for(auto j = _track_map->get(i)->begin_cluster_keys(); j!=_track_map->get(i)->end_cluster_keys(); j++){
      TrkrCluster *cl = _cluster_map->findCluster(*j);
      cout << cl->getX() << "," << cl->getY() << "," << cl->getZ() << endl;
    }

    }*/


  /*vector<vector<TrkrCluster*>> seedslist;
  vector<vector<double>> catalog;
  for(int r=0;r<60;r++){
    if(nlayer[r]>0 || r<50){
      printf("layer %d has radius %f\n",r,_radii_all[r]);
      }
  }
  PHTimer* t_fill2 = new PHTimer("t_fill2");
  t_fill2->stop();
  t_fill2->restart();
  PHTimer* t_assembly = new PHTimer("t_assembly");
  t_assembly->stop();
  t_assembly->restart();
  for(int r=2; r<52; r++){
    if(nlayer[r-2]>0&&nlayer[r-1]>0&&nlayer[r]>0&&nlayer[r+1]>0&&nlayer[r+2]>0)
      {
	std::vector<pointKey> layer_r_clusters;
	rtree.query(bgi::intersects(box(point(0,-3,r-0.5),point(2*M_PI,3,r+0.5))),std::back_inserter(layer_r_clusters));
	for(std::vector<pointKey>::iterator it_r=layer_r_clusters.begin(); it_r!=layer_r_clusters.end();it_r++)
	  {
	    std::vector<pointKey> in_clusters, out_clusters;
	    double this_phi = it_r->first.get<0>();
	    double this_eta = it_r->first.get<1>();
	    wquery(rtree,this_phi-phisr,this_eta-etasr,r-2.5,this_phi+phisr,this_eta+etasr,r-1.5,in_clusters);
	    wquery(rtree,this_phi-phisr,this_eta-etasr,r+1.5,this_phi+phisr,this_eta+etasr,r+2.5,out_clusters);
	    for(std::vector<pointKey>::iterator it_in=in_clusters.begin(); it_in!=in_clusters.end();it_in++)
	      {
		for(std::vector<pointKey>::iterator it_out=out_clusters.begin(); it_out!=out_clusters.end();it_out++)
		  {
		    double in_phi = it_in->first.get<0>();
		    double in_eta = it_in->first.get<1>();
		    double out_phi = it_out->first.get<0>();
		    double out_eta = it_out->first.get<1>();
		    if(abs(phidiff(this_phi,in_phi)-phidiff(out_phi,this_phi))>phist || abs(this_eta-in_eta-(out_eta-this_eta))>etast) {continue;}
		    std::vector<pointKey> xin_clusters, xout_clusters;
		    double exp_xin_phi = phiadd(this_phi,(_radii_all[r]-_radii_all[r-1])/(_radii_all[r]-_radii_all[r-2])*phidiff(in_phi,this_phi));
		    double exp_xout_phi = phiadd(this_phi,(_radii_all[r+1]-_radii_all[r])/(_radii_all[r+2]-_radii_all[r])*phidiff(out_phi,this_phi));
		    double exp_xin_eta = this_eta+(_radii_all[r]-_radii_all[r-1])/(_radii_all[r]-_radii_all[r-2])*(in_eta-this_eta);
		    double exp_xout_eta = this_eta+(_radii_all[r+1]-_radii_all[r])/(_radii_all[r+2]-_radii_all[r])*(out_eta-this_eta);
		    wquery(rtree,exp_xin_phi-phixt,exp_xin_eta-etaxt,r-1.5,exp_xin_phi+phixt,exp_xin_eta+etaxt,r-0.5,xin_clusters);
		    wquery(rtree,exp_xout_phi-phixt,exp_xout_eta-etaxt,r+0.5,exp_xout_phi+phixt,exp_xout_eta+etaxt,r+1.5,xout_clusters);
		    if(xin_clusters.empty() || xout_clusters.empty()) {continue;}
		    double xinrecord = 100.0, xoutrecord = 100.0;
		    pointKey *xinkey = &*xin_clusters.begin();
		    pointKey *xoutkey = &*xout_clusters.begin();
		    for(std::vector<pointKey>::iterator it=xin_clusters.begin(); it!=xin_clusters.end();it++){
		      double dist = abs(phidiff(it->first.get<0>(),exp_xin_phi))+abs(it->first.get<1>()-exp_xin_eta);
		      if(dist < xinrecord){
			*xinkey = *it;
			xinrecord = dist;
		      }
		    }
		    for(std::vector<pointKey>::iterator it=xout_clusters.begin(); it!=xout_clusters.end();it++){
		      double dist = abs(phidiff(it->first.get<0>(),exp_xout_phi))+abs(it->first.get<1>()-exp_xout_eta);
		      if(dist < xoutrecord){
			*xoutkey = *it;
			xoutrecord = dist;
		      }
		    }
		    vector<TrkrCluster*> mynewseed;
		    mynewseed.push_back(_cluster_map->findCluster(it_r->second));
		    mynewseed.push_back(_cluster_map->findCluster(it_in->second));
		    mynewseed.push_back(_cluster_map->findCluster(it_out->second));
		    mynewseed.push_back(_cluster_map->findCluster(xinkey->second));
		    mynewseed.push_back(_cluster_map->findCluster(xoutkey->second));
		    //vector<TrkrCluster*>* mynewseed = {new TrkrCluster(_cluster_map->findCluster(it_r->second))};
		    //vector<TrkrCluster*>* mynewseed = *{_cluster_map->findCluster(it_r->second), _cluster_map->findCluster(it_in->second), _cluster_map->findCluster(it_out->second), _cluster_map->findCluster(xinkey->second), _cluster_map->findCluster(xoutkey->second)};
		    sort((mynewseed).begin(),(mynewseed).end(),comparing);
		    seedslist.push_back(mynewseed);
		  }
	      }
	  }
      }
  }
  t_assembly->stop();
  cout << "This assembly time is " << t_assembly->get_accumulated_time()/1000 << " s" << endl;
  bool allmerged = false;
  cout << seedslist.size() << endl;
  //assert (seedslist.size() >= 2);
  while(!allmerged){
    cout << "not all merged" << endl;
    allmerged = true;
    if(seedslist.size() < 2) break;
    for(unsigned int i=0; i<seedslist.size()-1; i++){
      //allmerged = true;
      bool partiallymerged = false;
      while(!partiallymerged){
	partiallymerged = true;
	unsigned int j=i+1;
	while(j<seedslist.size()){
	  //partiallymerged = true;
	  if(std::includes(seedslist.at(i).begin(),seedslist.at(i).end(),seedslist.at(j).begin(),seedslist.at(j).end(),comparing)){
	  //if(issuperset(seedslist.at(i),seedslist.at(j))){
	  //if(equal(seedslist.at(i).begin(),seedslist.at(i).end(),yunion(seedslist.at(i),seedslist.at(j)).begin(),comparing)){
	    cout << "first condition invoked" << endl;
	    allmerged = false;
	    partiallymerged = false;
	    seedslist.erase(seedslist.begin()+j);
	  }
	  if(std::includes(seedslist.at(j).begin(),seedslist.at(j).end(),seedslist.at(i).begin(),seedslist.at(i).end(),comparing)){
	//else if(issuperset(seedslist.at(j),seedslist.at(i))){
	//else if(equal(seedslist.at(j).begin(),seedslist.at(j).end(),yunion(seedslist.at(i),seedslist.at(j)).begin(),comparing)){
	    cout << "second condition invoked" << endl;
	    allmerged = false;
	    partiallymerged = true;
	    seedslist.erase(seedslist.begin()+i);
	    j=seedslist.size();
	  }
	  else if(largeintersection((seedslist.at(i)),(seedslist.at(j)))){
	    allmerged = false;
	    partiallymerged = false;
	    changetounion((seedslist.at(i)),(seedslist.at(j)));
	    //delete seedslist.at(j);
	    seedslist.erase(seedslist.begin()+j);
	  }
	  else j++;
	}
      }
    }
  }
  cout << "seedslist has size " << seedslist.size() << endl;


  TH1D *phisrs = new TH1D("phisrs","phisrs",100,0,0.1);
  TH1D *etasrs = new TH1D("etasrs","etasrs",100,0,0.015);
  TH1D *phists = new TH1D("phisrs","phists",100,0,0.04);
  TH1D *etasts = new TH1D("etasts","etasts",100,0,0.01);
  TH1D *phixts = new TH1D("phixts","phixts",100,0,0.03);
  TH1D *etaxts = new TH1D("etaxts","etaxts",100,0,0.008);
  TH1D *times = new TH1D("times","times",100,0,200);
  double currenttime = 0;
  PHTimer* looptimer = new PHTimer("looptimer");




  for(unsigned int i=0; i<seedslist.size(); i++){
    looptimer->stop();
    if(i>0){times->Fill(looptimer->get_accumulated_time()-currenttime);currenttime=looptimer->get_accumulated_time();}
    looptimer->restart();
    clusterpoints = seedslist.at(i);
    cout << "Working on a seed of size " << clusterpoints.size() << endl;
    cout << "{";
    for(unsigned int j=0; j<clusterpoints.size(); j++){
      printf("{%d, %f, %lX}", TrkrDefs::getLayer(clusterpoints.at(j)->getClusKey()),atan2(clusterpoints.at(j)->getY(),clusterpoints.at(j)->getX()),clusterpoints.at(j)->getClusKey());
      //cout << "{" << TrkrDefs::getLayer(clusterpoints.at(j)->getClusKey()) << ", " << atan2(clusterpoints.at(j)->getY(),clusterpoints.at(j)->getX())<< "}";
      //cout << "{" << clusterpoints.at(j)->getX() << "," << clusterpoints.at(j)->getY() << "," << clusterpoints.at(j)->getZ() << "},";
    }
    cout << "}" << endl;
    ROOT::Math::Minimizer* min =
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
    
    // set tolerance , etc...
    min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    min->SetMaxIterations(10000);  // for GSL
    min->SetTolerance(0.001);
    min->SetPrintLevel(1);
    
    // create funciton wrapper for minmizer
    ROOT::Math::Functor f(&RosenBrock,5+clusterpoints.size());
    
    min->SetFunction(f);
    double Cx = clusterpoints.at(0)->getX();
    double Cy = clusterpoints.at(0)->getY();
    double Cz = clusterpoints.at(0)->getZ();
    min->SetVariable(0, "xc", Cx, 0.01);
    min->SetVariable(1, "yc", Cy, 0.01);
    min->SetVariable(2, "phi", atan2(Cy,Cx), 0.01);
    min->SetVariable(3, "cotheta", atan(Cz/sqrt(Cx*Cx+Cy*Cy)), 0.01);
    min->SetVariable(4, "kappa", 0.000001, 0.01);
    for(unsigned int j = 0; j < clusterpoints.size(); j++){
      min->SetVariable(5+j, Form("zhat%d", j),clusterpoints.at(j)->getZ(),0.01);
    }
    min->SetVariable(5+clusterpoints.size(),"zvertex",_vertex->get_z(),0.01);
    //cout <<  << ", " << 
    min->Minimize();
    if(min->Status()>1 || std::isnan(min->MinValue())) continue;
    
    for(unsigned int j = 2; j<clusterpoints.size()-2; j++){
      TVector3 vec0(clusterpoints.at(j-2)->getX()-_vertex->get_x(), clusterpoints.at(j-2)->getY()-_vertex->get_y(), clusterpoints.at(j-2)->getZ()-_vertex->get_z());
      TVector3 vec1(clusterpoints.at(j-1)->getX()-_vertex->get_x(), clusterpoints.at(j-1)->getY()-_vertex->get_y(), clusterpoints.at(j-1)->getZ()-_vertex->get_z());
      TVector3 vec2(clusterpoints.at(j)->getX()-_vertex->get_x(), clusterpoints.at(j)->getY()-_vertex->get_y(), clusterpoints.at(j)->getZ()-_vertex->get_z());
      TVector3 vec3(clusterpoints.at(j+1)->getX()-_vertex->get_x(), clusterpoints.at(j+1)->getY()-_vertex->get_y(), clusterpoints.at(j+1)->getZ()-_vertex->get_z());
      TVector3 vec4(clusterpoints.at(j+2)->getX()-_vertex->get_x(), clusterpoints.at(j+2)->getY()-_vertex->get_y(), clusterpoints.at(j+2)->getZ()-_vertex->get_z());

      phisrs->Fill(abs(phidiff(vec2.Phi(),vec0.Phi())));
      phisrs->Fill(abs(phidiff(vec2.Phi(),vec4.Phi())));
      etasrs->Fill(abs(vec2.Eta()-vec0.Eta()));
      etasrs->Fill(abs(vec2.Eta()-vec4.Eta()));
      phists->Fill(abs(phidiff(vec2.Phi(),vec0.Phi())-phidiff(vec4.Phi(),vec2.Phi())));
      etasts->Fill(abs((vec2.Eta()-vec0.Eta())-(vec4.Eta()-vec2.Eta())));


      double exp_xin_phi = phiadd(vec2.Phi(),(_radii_all[j]-_radii_all[j-1])/(_radii_all[j]-_radii_all[j-2])*phidiff(vec0.Phi(),vec2.Phi()));
      double exp_xout_phi = phiadd(vec2.Phi(),(_radii_all[j+1]-_radii_all[j])/(_radii_all[j+2]-_radii_all[j])*phidiff(vec4.Phi(),vec2.Phi()));
      double exp_xin_eta = vec2.Eta()+(_radii_all[j]-_radii_all[j-1])/(_radii_all[j]-_radii_all[j-2])*(vec0.Eta()-vec2.Eta());
      double exp_xout_eta = vec2.Eta()+(_radii_all[j+1]-_radii_all[j])/(_radii_all[j+2]-_radii_all[j])*(vec4.Eta()-vec2.Eta());
      phixts->Fill(abs(phidiff(vec1.Phi(),exp_xin_phi)));
      phixts->Fill(abs(phidiff(vec3.Phi(),exp_xout_phi)));
      etaxts->Fill(abs(vec1.Eta()-exp_xin_eta));
      etaxts->Fill(abs(vec3.Eta()-exp_xout_eta));

    }



    SvtxTrack_v1 track;
    track.set_id(numberofseeds);
    for(unsigned int j = 0; j<clusterpoints.size();j++) track.insert_cluster_key(clusterpoints.at(j)->getClusKey());
    track.set_chisq(min->MinValue());
    track.set_ndf(2*clusterpoints.size()-5);
    short int helicity = 1;
    const double *pars = min->X();
    double Qxc = pars[0];
    double Qyc = pars[1];
    double Qphi = pars[2];
    double Qcotheta = pars[3]; 
    double Qkappa = pars[4];
    if(Cz*Qcotheta<0){
      Qphi *=-1;
      Qcotheta *=-1;
      Qkappa *=-1;
    }
    vector<double> insertion{Qxc,Qyc,Qphi,Qcotheta,Qkappa,Cz};
    catalog.push_back(insertion);

    double CVxcxc = min->CovMatrix(0,0);
    double CVxcyc = min->CovMatrix(0,1);
    //double CVxczc = min->CovMatrix(0,5);
    double CVycyc = min->CovMatrix(1,1);
    //double CVyczc = min->CovMatrix(1,5);
    //double CVzczc = min->CovMatrix(5,5);
    double CVphiphi = min->CovMatrix(2,2);
    double CVphicotheta = min->CovMatrix(2,3);
    double CVphikappa = min->CovMatrix(2,4);
    double CVcothetacotheta = min->CovMatrix(3,3);
    double CVcothetakappa = min->CovMatrix(3,4);
    double CVkappakappa = min->CovMatrix(4,4);
    const double BQ = 0.01*1.4*0.299792458;//1.4 / 333.6;
    if (Qcotheta*Qkappa<0) helicity -= 2;
    track.set_charge(helicity);
    track.set_x(Qxc);
    track.set_y(Qyc);
    track.set_z(Cz);
    track.set_px(BQ/Qkappa*cos(Qphi));
    track.set_py(BQ/Qkappa*sin(Qphi));
    track.set_pz(BQ/Qkappa*tan(Qcotheta));
    double BQ2 = BQ*BQ;
    double CP = cos(Qphi);
    double SP = sin(Qphi);
    double K2 = Qkappa*Qkappa;
    double K3 = K2*Qkappa;
    double K4 = K3*Qkappa;
    double TT = tan(Qcotheta);
    double ST = 1/(cos(Qcotheta)*cos(Qcotheta));
    track.set_error(0,0,CVxcxc);
    track.set_error(1,1,CVycyc);
    track.set_error(2,2,0);
    track.set_error(0,1,CVxcyc);
    track.set_error(1,0,CVxcyc);
    track.set_error(0,2,0);
    track.set_error(2,0,0);
    track.set_error(1,2,0);
    track.set_error(2,1,0);
    track.set_error(3,3,BQ2*(CP*CP/K4*CVkappakappa+2*CP*SP/K3*CVphikappa+SP*SP/K2*CVphiphi));
    track.set_error(4,4,BQ2*(CP*CP/K2*CVphiphi-2*CP*SP/K3*CVphikappa+SP*SP/K4*CVkappakappa));
    track.set_error(5,5,BQ2*(ST*ST/K2*CVcothetacotheta-2*ST*TT/K3*CVcothetakappa+TT*TT/K4*CVkappakappa));
    track.set_error(3,4,BQ2*(-CP*CP/K3*CVphikappa+CP*SP/K4*CVkappakappa-CP*SP/K2*CVphiphi+SP*SP/K3*CVphikappa));
    track.set_error(4,3,BQ2*(-CP*CP/K3*CVphikappa+CP*SP/K4*CVkappakappa-CP*SP/K2*CVphiphi+SP*SP/K3*CVphikappa));
    track.set_error(3,5,BQ2*(-CP*ST/K3*CVcothetakappa-SP*ST/K2*CVphicotheta+CP*TT/K4*CVkappakappa+SP*TT/K3*CVphikappa));
    track.set_error(5,3,BQ2*(-CP*ST/K3*CVcothetakappa-SP*ST/K2*CVphicotheta+CP*TT/K4*CVkappakappa+SP*TT/K3*CVphikappa));
    track.set_error(4,5,BQ2*(-SP*ST/K3*CVcothetakappa+CP*ST/K2*CVphicotheta+SP*TT/K4*CVkappakappa-CP*TT/K3*CVphikappa));
    track.set_error(5,4,BQ2*(-SP*ST/K3*CVcothetakappa+CP*ST/K2*CVphicotheta+SP*TT/K4*CVkappakappa-CP*TT/K3*CVphikappa));
    printf("Labeled[Show[ListPointPlot3D[Import[\"/Users/robert/Desktop/test.csv\"]], generateplot[%f,%f,%f,%f,%f,%f,-20, 5,Lighter[Red]]],{\"\\[Chi]^2/ndf=%f\\npx=%f\\[PlusMinus]%f\\npy=%f\\[PlusMinus]%f\\npz=%f\\[PlusMinus]%f\"},{Right}]\n",Qxc,Qyc,Qphi,Qcotheta,Qkappa,Cz,min->MinValue()/(2*clusterpoints.size()-5),BQ/Qkappa*cos(Qphi),sqrt(BQ2*(CP*CP/K4*CVkappakappa+2*CP*SP/K3*CVphikappa+SP*SP/K2*CVphiphi)),BQ/Qkappa*sin(Qphi),sqrt(BQ2*(CP*CP/K2*CVphiphi-2*CP*SP/K3*CVphikappa+SP*SP/K4*CVkappakappa)),BQ/Qkappa*tan(Qcotheta),sqrt(BQ2*(ST*ST/K2*CVcothetacotheta-2*ST*TT/K3*CVcothetakappa+TT*TT/K4*CVkappakappa)));
    for(unsigned int par1=0; par1<3; par1++)
      for(unsigned int par2 = 0; par2<3; par2++)
	{
	  track.set_error(par1,par2+3,0);
	  track.set_error(par1+3,par2,0);
	}
	_track_map->insert(&track);
    numberofseeds +=1;
  }
  cout << "Number of seeds: " << numberofseeds << endl;
  for(auto i=catalog.begin();i!=catalog.end();i++)
  {
    //cout << *i[0] << "," << *i[1] << "," << *i[2] << "," << *i[3] << "," << *i[4] << "," << *i[5] << endl;
    cout << i->at(0) << "," << i->at(1) << "," << i->at(2) << "," << i->at(3) << "," << i->at(4) << "," << i->at(5) << endl;
  }
  t_fill2->stop();
  cout << "seeding time is " << t_fill2->get_accumulated_time()/1000 << " s" << endl;
  phisrs->Write("phisrs");
  etasrs->Write("etasrs");
  phists->Write("phists");
  etasts->Write("etasts");
  phixts->Write("phixts");
  etaxts->Write("etaxts");
  times->Write("times");*/

  return Fun4AllReturnCodes::EVENT_OK;
  }


int PHRTreeSeeding::Setup(PHCompositeNode* topNode)
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


