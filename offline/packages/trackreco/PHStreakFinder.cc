/*!
 *  \file PHStreakFinder.cc
 
 *  \author Christof Roland 
 */

//begin

//#ifndef PHCOSMICSFILTER_H
//#define PHCOSMICSFILTER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>

// Helix Hough + Eigen includes (hidden from rootcint)
#if !defined(__CINT__) || defined(__CLING__)

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

#include "PHStreakFinder.h"

#include "AssocInfoContainer.h"                         // for AssocInfoCont...

// trackbase_historic includes
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>  // for getLayer, clu...
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase/LaserCluster.h>
#include <trackbase/LaserClusterContainer.h>
#include <trackbase/LaserClusterContainerv1.h>
#include <trackbase/LaserClusterv1.h>

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



//ROOT includes for debugging
#include <TFile.h>
#include <TNtuple.h>
#include <TAxis.h>
#include <TGraph.h>
//BOOST for combi seeding
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

using point = bg::model::point<float, 3, bg::cs::cartesian>;
//typedef bg::model::point<float, 3, bg::cs::cartesian> point;
using box = bg::model::box<point>;
//typedef bg::model::box<point> box;
using pointKey = std::pair<point, int>;
using stubKey = std::pair<point, float>;
//typedef std::pair<point, TrkrDefs::cluskey> pointKey;

using pointrtree = bgi::rtree<pointKey, bgi::quadratic<16>>;
using stubrtree = bgi::rtree<pointKey, bgi::quadratic<16>>;

// standard includes
#include <TH1.h>
//#include <stdio.h>      /* printf */
//#include <cmath.h>       /* copysign */
#include <algorithm>
#include <climits>                                     // for UINT_MAX
#include <cmath>
#include <iostream>
#include <memory>
#include <set>                                          // for set
#include <tuple>
#include <utility>   
                                   // for pair, make_pair

#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp


using namespace std;
//using namespace ROOT::Minuit2;
//namespace bg = boost::geometry;
//namespace bgi = boost::geometry::index;

//typedef uint64_t cluskey;

//vector<LaserCluster*> clusterpoints;
using myvec = vector<float>;

PHStreakFinder::PHStreakFinder(const std::string& name, const string& filename)
  : SubsysReco(name)
, _filename(filename)
{
  
}

int PHStreakFinder::GetNodes(PHCompositeNode* topNode)
{
  tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tGeometry)
  {
    std::cout << PHWHERE << "No acts reco geometry, bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _cluster_map = findNode::getClass<LaserClusterContainerv1>(topNode, "LASER_CLUSTER");
  if (!_cluster_map)
  {
    std::cout << PHWHERE << "No cluster container on node tree, bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHStreakFinder::Init(PHCompositeNode* topNode)
{
  if(topNode==nullptr){
    std::cout << PHWHERE << "No topNode, bailing."
	      << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _nevent = 0;
    _tfile = new TFile(_filename.c_str(), "RECREATE");
  _ntp_ev = new TNtuple("ntp_ev", "cos event info","ev:ev2");
  _ntp_clu = new TNtuple("ntp_clu", "cos event info","ev:x:y:z");
  _ntp_clutrk = new TNtuple("ntp_clutrk", "cos event info","ev:ntrk:x:y:z:nclus:nclus0:nclus1:r:phi:dedx:intxz:intyz:intrz:slxz:slyz:slrz:mresx:mresy:zmin:zmax");
  _ntp_stub = new TNtuple("ntp_stub", "cos stub info","ev:intxz:slxz:intyz:slyz:tanxz:tanyz:num");
  _ntp_trk = new TNtuple("ntp_trk", "cos stub info","ev:ntrk:nclus:nclus0:nclus1:r:phi:dedx:intxz:intyz:intrz:slxz:slyz:slrz:mresx:mresy:zmin:zmax");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHStreakFinder::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }
  ret = createNodes(topNode);


  return ret;
}
int PHStreakFinder::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHStreakFinder::createNodes");
  }

  PHNodeIterator dstIter(dstNode);
  PHCompositeNode* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_seedContainer = findNode::getClass<TrackSeedContainer>(topNode, m_trackMapName);
  if (!m_seedContainer)
  {
    m_seedContainer = new TrackSeedContainer_v1;
    PHIODataNode<PHObject>* trackNode =
        new PHIODataNode<PHObject>(m_seedContainer, m_trackMapName, "PHObject");
    svtxNode->addNode(trackNode);
  }
  if (m_seedContainer){
    cout << "SEED CONTAINER CREATED" << endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

double PHStreakFinder::phiadd(double phi1, double phi2){
  double s = phi1+phi2;
  if(s>2*M_PI) {return s-2*M_PI;}
  else if(s<0) {return s+2*M_PI;}
  else {return s;}
}


double PHStreakFinder::phidiff(double phi1, double phi2){
  double d = phi1-phi2;
  if(d>M_PI) {return d-2*M_PI;}
  else if(d<-M_PI) {return d+2*M_PI;}
  else {return d;}
}

void PHStreakFinder::get_stub(const pointrtree &search_rtree, float pointx, float pointy, float pointz, int &count, double &xzslope, double &xzintercept, double &yzslope, double &yzintercept){ //NOLINT
  float m1_dx = 2;
  float m1_dy = 2;
  float m1_dz = 8;
  vector<pointKey> boxclusters;
  search_rtree.query(bgi::intersects(box(point(pointx-m1_dx,pointy-m1_dy,pointz-m1_dz),
				  point(pointx+m1_dx,pointy+m1_dy,pointz+m1_dz))),std::back_inserter(boxclusters));
  
  int nbox = boxclusters.size();
  
  int nhit = 0;
  double xsum = 0;
  double x2sum = 0;
  double ysum = 0;
  double y2sum = 0;
  double z2sum = 0;

  double zsum = 0;
  double xzsum = 0;
  double yzsum = 0;

  if(nbox>=2){
    for(vector<pointKey>::iterator pbox = boxclusters.begin();pbox!=boxclusters.end();++pbox){
      float boxx = pbox->first.get<0>();
      float boxy = pbox->first.get<1>();
      float boxz = pbox->first.get<2>();
      //      std::cout << "fprobe: " << count << " x:"  << boxx << " y: " << boxy << std::endl;
      nhit++;
      float x = boxx;
      float y = boxy;
      float z = boxz;

      xsum = xsum + x;            // calculate sigma(xi)
      ysum = ysum + y;            // calculate sigma(xi)
      zsum = zsum + z;            // calculate sigma(yi)

      x2sum = x2sum + x*x;  // calculate sigma(x^2i)
      y2sum = y2sum + y*y;  // calculate sigma(x^2i)
      z2sum = z2sum + z*z;  // calculate sigma(x^2i)

      xzsum = xzsum + x * z;      // calculate sigma(xi*yi)
      yzsum = yzsum + y * z;      // calculate sigma(xi*yi)
    }
  }
  
  const double denominator = ((z2sum * nhit) - (zsum*zsum));
  //  const double denominator_yz = ((y2sum * nhit) - (ysum*ysum));

  xzslope = (xzsum * nhit - zsum * xsum) / denominator;   // calculate slope
  xzintercept = (z2sum * xsum - zsum * xzsum) / denominator;  // calculate intercept

  yzslope = (yzsum * nhit - zsum * ysum) / denominator;   // calculate slope
  yzintercept = (z2sum * ysum - zsum * yzsum) / denominator;  // calculate intercept

  count = nhit;
}

int PHStreakFinder::process_event(PHCompositeNode* topNode)
{
  if(topNode==nullptr){
    std::cout << PHWHERE << "No topNode, bailing."
	      << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
    
  _nevent++;
  _ntp_ev->Fill(_nevent,_nevent);
  int ntrack = 0;
  std::vector<TrackSeed_v2> clean_chains;
  //Fill rtree
  bgi::rtree<pointKey, bgi::quadratic<16> > rtree;

  auto range = _cluster_map->getClusters();
  for( auto clusIter = range.first; clusIter != range.second; ++clusIter ){
    LaserCluster *cluster = clusIter->second;
    _ntp_clu->Fill(_nevent,cluster->getX() , cluster->getY(), cluster->getZ());
    rtree.insert(std::make_pair(point(cluster->getX() , cluster->getY(), cluster->getZ()),std::distance(range.first,clusIter) ));
      
  }

  //Get all clusters from rtree, fit stubs around clusters, fill stub tree
  // search horizontal strips in x - y search grid

  for(float ix = -80;ix<80;ix+=4){
    for(float iy = -80;iy<80;iy+=4){ 
      float ir = sqrt(ix*ix+iy*iy);
      if(ir<30||ir>78){
	continue;
      }
      vector<pointKey> allclusters;
      
      rtree.query(bgi::intersects(box(point(ix-2,iy-2,-110),point(ix+2,iy+2,110))),std::back_inserter(allclusters));
      //if(Verbosity()>0)
      if(allclusters.size()>20){
	//cout << " ev: " << _nevent << " ix = " << ix << " iy: " << iy << " number clus is " << allclusters.size() << endl;
      
	
	bgi::rtree<stubKey, bgi::quadratic<16> > xzrtree_stub;
	bgi::rtree<stubKey, bgi::quadratic<16> > yzrtree_stub;
	
	for(vector<pointKey>::iterator cluster = allclusters.begin();cluster!=allclusters.end();++cluster){
	  float pointx = cluster->first.get<0>();
	  float pointy = cluster->first.get<1>();
	  float pointz = cluster->first.get<2>();
	  int fcount = 0;
	  double xzslope = 0;
	  double xzintercept = 0;
	  double yzslope = 0;
	  double yzintercept = 0;
	  //calc slope and intersect from 5 cluster stubs
	  get_stub(rtree, pointx, pointy, pointz, fcount, xzslope, xzintercept, yzslope, yzintercept);
	  //int num = 0;
	  if(finite(xzslope)&&finite(yzslope)){
	    xzrtree_stub.insert(std::make_pair(point(xzintercept , xzslope, yzintercept), yzslope));
	    //  cout << " inserting: " << "xzint: " << xzintercept << "xzsl: " << xzslope << " yzint: " << yzintercept << " yzsl: " << yzslope <<endl;
	    float tana_xz = atan(xzslope);
	    float tana_yz = atan(yzslope);
	    //ev:intxz:slxz:intyz:slyz:tanxz:tanyz:xo
	    if(_write_ntp)_ntp_stub->Fill(_nevent,xzintercept,xzslope,yzintercept,yzslope,tana_xz,tana_yz,fcount);
	    
	  }	  
	}
	//find good stub candidates
	map <int,pair<float,float>> xzouttrkmap;
	map <int,myvec> xztrkmap;
	//	int nout_xz = 0;
	float intmax = 80;
	float intmin = -80;
	float slmax = 0.2;
	float slmin = -0.2;
	float int_width = (intmax - intmin)/20;
	float sl_width = (slmax - slmin)/8;
	cout << "number xzstub is " << xzrtree_stub.size() << endl;
	
	std::vector<stubKey> xzallstubs;
	xzrtree_stub.query(bgi::intersects(box(point(intmin,slmin,intmin),point(intmax,slmax,intmax))),std::back_inserter(xzallstubs));
	if(xzallstubs.size()==0){cout<<"nostubs" <<endl; break;}


	for(vector<stubKey>::iterator stub = xzallstubs.begin();stub!=xzallstubs.end();++stub){
	  float pxzint = stub->first.get<0>();
	  float pxzsl = stub->first.get<1>();
	  float pyzint = stub->first.get<2>();
	  float pyzsl =  stub->second;
	  //cout << "xzint: " << pxzint << "xzsl: " << pxzsl << " yzint: " << pyzint << " yzsl: " << pyzsl <<endl;
	  if(abs(pyzsl)<0.2){
	    vector<stubKey> trkcand;
	    xzrtree_stub.query(bgi::intersects(box(point(pxzint-int_width,pxzsl-sl_width,pyzint-int_width),point(pxzint+int_width,pxzsl+sl_width,pyzint+int_width))),std::back_inserter(trkcand));
      
	    int ntrk = trkcand.size();
	    int count = 0;
	    float xzintsum = 0;
	    float xzslsum = 0;
	    float yzintsum = 0;
	    float yzslsum = 0;
	    if(ntrk>=5){
	      for(vector<stubKey>::iterator ptrk = trkcand.begin();ptrk!=trkcand.end();++ptrk){
		float trkxzint = ptrk->first.get<0>();
		float trkxzsl = ptrk->first.get<1>();
		float trkyzint = ptrk->first.get<2>();
		float trkyzsl = ptrk->second;
		//	  if(Verbosity()>0)
		/*cout<< " ev: " << _nevent
		    << " stub " << ntrk
		    << " xzint: " << trkxzint 
		    << " xzsl: "  << trkxzsl 
		    << " yzint: " << trkyzint 
		    << " yzsl: "  << trkyzsl 
		    << endl;
		*/
		xzintsum = xzintsum + trkxzint;            // calculate sigma(xi)
		xzslsum = xzslsum + trkxzsl;            // calculate sigma(yi)
		yzintsum = yzintsum + trkyzint;            // calculate sigma(xi)
		yzslsum = yzslsum + trkyzsl;            // calculate sigma(yi)
		count++;
	      }
	      float mxzint = (xzintsum/count);
	      float mxzsl  = (xzslsum/count);
	      float myzint = (yzintsum/count);
	      float myzsl  = (yzslsum/count);
	      myvec outarr = {mxzint,mxzsl,myzint,myzsl} ;
	      /*	      outarr.push_back(mxzint);
	      outarr.push_back(mxzsl);
	      outarr.push_back(myzint);
	      outarr.push_back(myzsl);
	      */
	      xztrkmap[ntrk] = outarr;
	      //.insert(std::pair<int,vector<float>(ntrk,outarr)); 
	      
	      if(Verbosity()>0) {
	      cout<< " sxz tub in box " << ntrk
		  << " xzint: " << pxzint 
		  << " xzsl: " << pxzsl 
		  << " yzint: " << pyzint 
		  << " yzsl: " << pyzsl 
		  << " intwi: " << int_width 
		  << " slwi: " << sl_width
		  << " mxzint " << mxzint
		  << " mxzsl " << mxzsl
		  << " myzint " << myzint
		  << " myzsl " << myzsl
		  << endl;
	      }
	    }
	  }
	}
	float meanxzint = 0;
	float meanxzsl = 0;
	float meanyzint = 0;
	float meanyzsl = 0;
//	int meancnt = 0;
	for(auto ptrk = xztrkmap.begin();ptrk!=xztrkmap.end();++ptrk){
	  int nstub = ptrk->first;
	  myvec value = ptrk->second;
//	  meancnt++;
	  meanxzint = value[0];
	  meanxzsl  = value[1];
	  meanyzint = value[2];
	  meanyzsl  = value[3];
	  //	  for (const auto& [nstub, value[4]] : xztrkmap)
	   if(Verbosity()>0) {
	     std::cout <<" nstub: " << nstub
		    << "xzint " << value[0] 
		    << "xzsl  " << value[1] 
		    << "yzint " << value[2] 
		    << "yzsl " << value[3] 
		    << std::endl;
	   }
	}
	//collect hits and fill ntuple
	
	vector<pointKey> trkclusters;
	vector<pointKey> lineclusters;
	rtree.query(bgi::intersects(box(point(-80,-80,-110),point(80,80,110))),std::back_inserter(lineclusters));
	// if(Verbosity()>0) 
	 if(Verbosity()>0) {
	   cout << "number line clus is " << lineclusters.size() << endl;
	 }
	float meanresx = 0;
	float meanresy = 0;
	  
	for(vector<pointKey>::iterator clustertrk = lineclusters.begin();clustertrk!=lineclusters.end();++clustertrk){
      
	  float ptx = clustertrk->first.get<0>();
	  float pty = clustertrk->first.get<1>();
	  float ptz = clustertrk->first.get<2>();
	  
	  float refx = meanxzint + ptz*meanxzsl;
	  float refy = meanyzint + ptz*meanyzsl;
	  float resx =  abs(ptx-refx);
	  float resy =  abs(pty-refy);
	  //      if(Verbosity()>0)
	  if((resx<2)&&(resy<2)){
	    meanresx += resx;
	    meanresy += resy;
	    if(Verbosity()>0) { cout << " x: " << ptx << " | " << refx << " y: " << pty << " | " << refy << " resx " << resx << " resy:" << resy << endl;}
	    trkclusters.push_back(*clustertrk);
	  }
	}
	meanresx/= trkclusters.size();
	meanresy/= trkclusters.size();
	int nclus = 0;
	int nclus0 = 0;
	int nclus1 = 0;
	float zmin = 10000;
	float zmax = -10000;
	std::vector<float> dedxlist;
	
	//cout << "number trk clus is " << trkclusters.size() << endl;
	if(trkclusters.size()>=_min_nclusters){
	  //cout << "setting keep true: " << trkclusters.size() << " > " << _min_nclusters << endl;
	  for (const auto& trkclu : trkclusters){
	    //      std::cout <<" yyz cand out ev: " << _nevent << '[' << xzkey << "] = " << yzvalue.first << " | " << xzvalue.second << std::endl;
	    // for(vector<pointKey>::iterator clustertrk = trkclusters.begin();clustertrk!=trkclusters.end();++clustertrk){
	    float tptx = trkclu.first.get<0>();
	    float tpty = trkclu.first.get<1>();
	    float tptz = trkclu.first.get<2>();
	    //      float ptx = clustertrk->first.get<0>();
	    //float pty = clustertrk->first.get<1>();
	    // float ptz = clustertrk->first.get<2>();
	    nclus++;
	    if(tptz<0){
	      nclus0++;
	    }
	    if(tptz>0){
	      nclus1++;
	    }
	    if(tptz<zmin){
	      zmin = tptz;
	    }
	    if(tptz>zmax){
	      zmax = tptz;
	    }
	    float adc = trkclu.second;
	    
	    float thick = 1.098;
	    float r  = sqrt(tptx*tptx+tpty*tpty);
	    if(r>30&&r<=40.5){
	      thick = 0.57;
	    }else if(r>40.5&&r<60){
	      thick = 1.02;
	    } 
	    dedxlist.push_back(adc/thick);
	    //calc dedx
	    //_ntp_clutrk->Fill(_nevent,tptx , tpty, tptz);
	  }
	  sort(dedxlist.begin(), dedxlist.end());
	  int trunc_min = 0;
	  int trunc_max = (int)dedxlist.size()*0.7;
	  float sumdedx = 0;
	  int ndedx = 0;
	  for(int j = trunc_min; j<=trunc_max;j++){
	    sumdedx+=dedxlist.at(j);
	    ndedx++;
	  }
	  sumdedx/=ndedx;
	  TVector2 vec(meanxzint,meanyzint);
	  // ev:nclus:nclus0:nclus1:r:phi:dedx:intxz:intyz:intrz:slxz:slyz:slrz
	  float r = sqrt(meanxzint*meanxzint+ meanyzint*meanyzint);
	  float xp = meanxzint+meanxzsl*100;
	  float xn = meanxzint+meanxzsl*-100;
	  float yp = meanyzint+meanyzsl*100;
	  float yn = meanyzint+meanyzsl*-100;
	  float rp = sqrt(xp*xp+yp*yp);
	  float rn = sqrt(xn*xn+yn*yn);
	  float phi = vec.Phi();
	  float dr = rp-rn;
	  float ftX[30] = {0};
	  int nft = 0;
	  ftX[nft++] = _nevent;
	  ftX[nft++] =ntrack;
	  ftX[nft++] =nclus;
	  ftX[nft++] =nclus0;
	  ftX[nft++] =nclus1;
	  ftX[nft++] =r;
	  ftX[nft++] =phi;
	  ftX[nft++] =sumdedx;
	  ftX[nft++] =meanxzint;
	  ftX[nft++] =meanyzint;
	  ftX[nft++] =r;
	  ftX[nft++] =meanxzsl;
	  ftX[nft++] =meanyzsl;
	  ftX[nft++] =dr/200;
	  ftX[nft++] =meanresx;
	  ftX[nft++] =meanresy;
	  ftX[nft++] =zmin;
	  ftX[nft++] =zmax;
	  _ntp_trk->Fill(ftX);

	  for (const auto& trkclu : trkclusters){
	    float tptx = trkclu.first.get<0>();
	    float tpty = trkclu.first.get<1>();
	    float tptz = trkclu.first.get<2>();
	    float fcX[30] = {0};
	    int nfc = 0;
	    fcX[nfc++] = _nevent;
	    fcX[nfc++] =ntrack;
	    fcX[nfc++] =tptx;
	    fcX[nfc++] =tpty;
	    fcX[nfc++] =tptz;
	    fcX[nfc++] =nclus;
	    fcX[nfc++] =nclus0;
	    fcX[nfc++] =nclus1;
	    fcX[nfc++] =r;
	    fcX[nfc++] =phi;
	    fcX[nfc++] =sumdedx;
	    fcX[nfc++] =meanxzint;
	    fcX[nfc++] =meanyzint;
	    fcX[nfc++] =r;
	    fcX[nfc++] =meanxzsl;
	    fcX[nfc++] =meanyzsl;
	    fcX[nfc++] =dr/200;
	    fcX[nfc++] =meanresx;
	    fcX[nfc++] =meanresy;
	    fcX[nfc++] =zmin;
	    fcX[nfc++] =zmax;
	    
	    _ntp_clutrk->Fill(fcX);
	  }
	  ntrack++;
	}
      }
    }
  }


 #ifdef KAKAKA

  
  //  int out = 0;
  // find clusters of slope intectept pairs
  map <int,pair<float,float>> xzouttrkmap;
  map <int,pair<float,float>> yzouttrkmap;
  map <int,pair<float,float>> xztrkmap;
  map <int,pair<float,float>> yztrkmap;
  int nout_xz = 0;
  int nout_yz = 0;
  float intmax = 80;
  float intmin = -80;
  float slmax = 0.5;
  float slmin = -0.5;
  float int_width = (intmax - intmin)/20;
  float sl_width = (slmax - slmin)/20;
  cout << "number xzstuc is " << xzrtree_stub.size() << endl;
  cout << "number yzstub is " << yzrtree_stub.size() << endl;
  while(xzrtree_stub.size()>10){
    std::vector<pointKey> xzallstubs;
    xzrtree_stub.query(bgi::intersects(box(point(intmin,slmin,-1),point(intmax,slmax,1))),std::back_inserter(xzallstubs));
    if(xzallstubs.size()==0){cout<<"nostubs" <<endl; break;}


    for(vector<pointKey>::iterator stub = xzallstubs.begin();stub!=xzallstubs.end();++stub){
      float pint = stub->first.get<0>();
      float psl = stub->first.get<1>();
      cout << "int: " << pint << "sl: " << psl << endl;
      vector<pointKey> trkcand;
      
      xzrtree_stub.query(bgi::intersects(box(point(pint-int_width,psl-sl_width,-1),point(pint+int_width,psl+sl_width,1))),std::back_inserter(trkcand));
      
      int ntrk = trkcand.size();
      int count = 0;
      float intsum = 0;
      float slsum = 0;
      if(ntrk>=5){
	for(vector<pointKey>::iterator ptrk = trkcand.begin();ptrk!=trkcand.end();++ptrk){
	  float trkint = ptrk->first.get<0>();
	  float trksl = ptrk->first.get<1>();
	  //	  if(Verbosity()>0)
	    cout<< " ev: " << _nevent
		<< " xz stub " << ntrk
		<< " int: " << trkint 
		<< " sl: "  << trksl 
		<< endl;
	  
	  intsum = intsum + trkint;            // calculate sigma(xi)
	  slsum = slsum + trksl;            // calculate sigma(yi)
	  count++;
	}
	float mint = (intsum/count);
	float msl  = (slsum/count);
	xztrkmap[ntrk] = std::make_pair(mint,msl); 
      
	//	if(Verbosity()>0) 
	cout<< " sxz tub in box " << ntrk
	    << " int: " << pint 
	    << " sl: " << psl 
	    << " intwi: " << int_width 
	    << " slwi: " << sl_width
	    << " mint " << mint
	    << " msl " << msl
	    << endl;
      }
    }
  while(yzrtree_stub.size()>10){
    std::vector<pointKey> yzallstubs;
    yzrtree_stub.query(bgi::intersects(box(point(intmin,slmin,-1),point(intmax,slmax,1))),std::back_inserter(yzallstubs));
    if(yzallstubs.size()==0){break;}

    for(vector<pointKey>::iterator stub = yzallstubs.begin();stub!=yzallstubs.end();++stub){
      float pint = stub->first.get<0>();
      float psl = stub->first.get<1>();
      vector<pointKey> trkcand;
      
      yzrtree_stub.query(bgi::intersects(box(point(pint-int_width,psl-sl_width,-1),point(pint+int_width,psl+sl_width,1))),std::back_inserter(trkcand));
      
      int ntrk = trkcand.size();
      int count = 0;
      float intsum = 0;
      float slsum = 0;
      if(ntrk>=5){
	for(vector<pointKey>::iterator ptrk = trkcand.begin();ptrk!=trkcand.end();++ptrk){
	  float trkint = ptrk->first.get<0>();
	  float trksl = ptrk->first.get<1>();
	  //	  if(Verbosity()>0)
	    cout<< " ev: " << _nevent
		<< " yz stub " << ntrk
		<< " int: " << trkint 
		<< " sl: "  << trksl 
		<< endl;
	  
	  intsum = intsum + trkint;            // calculate sigma(xi)
	  slsum = slsum + trksl;            // calculate sigma(yi)
	  count++;
	}
	float mint = (intsum/count);
	float msl  = (slsum/count);
	yztrkmap[ntrk] = std::make_pair(mint,msl); 
      
	//	if(Verbosity()>0) 
	cout<< " syz tub in box " << ntrk
	    << " int: " << pint 
	    << " sl: " << psl 
	    << " intwi: " << int_width 
	    << " slwi: " << sl_width
	    << " mint " << mint
	    << " msl " << msl
	    << endl;
      }
    }
    
    if( xztrkmap.size()>0){
      int size_before = xzrtree_stub.size();
      //      if(Verbosity()>0)
      cout << "mapend: " << xztrkmap.rbegin()->first << " | " << (xztrkmap.rbegin()->second).first << " | "  << xztrkmap.rbegin()->second.second << endl;
      // for (const auto& [key, value] : trkmap)
      //	std::cout <<" ev: " << _nevent << '[' << key << "] = " << value.first << " | " << value.second << std::endl;
      //remove stubs in box
      float rmint = xztrkmap.rbegin()->second.first;
      float rmsl = xztrkmap.rbegin()->second.second;
      vector<pointKey> rmcand;
      xzrtree_stub.query(bgi::intersects(box(point(rmint-int_width,rmsl-sl_width,-1),point(rmint+int_width,rmsl+sl_width,1))),std::back_inserter(rmcand));
      for(vector<pointKey>::iterator rmstub = rmcand.begin();rmstub!=rmcand.end();++rmstub){
	float rmpint = rmstub->first.get<0>();
	float rmpsl = rmstub->first.get<1>();
	if(Verbosity()>0) cout<< "    rm " <<  " int: " << rmpint 
			      << " sl: " << rmpsl 
			      << endl;
	
	//auto rmpoint = std::make_pair(point(rmpint , rmpsl, 0.0), 0);
	//rtree_stub.remove(rmpoint);
	//	pointKey match;
	xzrtree_stub.remove(*rmstub);
	/*	if (!rtree_stub.query(bgi::nearest(rmpoint, 1), &match))
		continue;
	*/
	/*
	  if (bg::distance(rmpoint, match.first) <= 0.5) {
	  rtree_stub.remove(match);
	  
	*/
      }
      int size_after = xzrtree_stub.size();
      if(size_before == size_after) {break;}
      xzouttrkmap[nout_xz++] = std::make_pair(rmint,rmsl);
      cout << " tree sixe after remove: " << xzrtree_stub.size() << endl;
    }else{
      break;
    }
    xztrkmap.clear();
  }
  if( yztrkmap.size()>0){
    int size_before = yzrtree_stub.size();
    //if(Verbosity()>0)
    cout << "mapend: " << yztrkmap.rbegin()->first << " | " << (yztrkmap.rbegin()->second).first << " | "  << yztrkmap.rbegin()->second.second << endl;
    // for (const auto& [key, value] : trkmap)
    //	std::cout <<" ev: " << _nevent << '[' << key << "] = " << value.first << " | " << value.second << std::endl;
    //remove stubs in box
    float rmint = yztrkmap.rbegin()->second.first;
    float rmsl = yztrkmap.rbegin()->second.second;
    vector<pointKey> rmcand;
    yzrtree_stub.query(bgi::intersects(box(point(rmint-int_width,rmsl-sl_width,-1),point(rmint+int_width,rmsl+sl_width,1))),std::back_inserter(rmcand));
    for(vector<pointKey>::iterator rmstub = rmcand.begin();rmstub!=rmcand.end();++rmstub){
      float rmpint = rmstub->first.get<0>();
      float rmpsl = rmstub->first.get<1>();
      if(Verbosity()>0) cout<< "    rm " <<  " int: " << rmpint 
			    << " sl: " << rmpsl 
			    << endl;
      
      //auto rmpoint = std::make_pair(point(rmpint , rmpsl, 0.0), 0);
      //rtree_stub.remove(rmpoint);
      //	pointKey match;
      yzrtree_stub.remove(*rmstub);
      /*	if (!rtree_stub.query(bgi::nearest(rmpoint, 1), &match))
		continue;
      */
      /*
	if (bg::distance(rmpoint, match.first) <= 0.5) {
	rtree_stub.remove(match);
	
      */
    }
    int size_after = yzrtree_stub.size();
    if(size_before == size_after) {break;}
    yzouttrkmap[nout_yz++] = std::make_pair(rmint,rmsl);
    cout << " tree sixe after remove: " << yzrtree_stub.size() << endl;
  }else{
    break;
  }
  yztrkmap.clear();
  }
  
  for (const auto& [key, value] : xzouttrkmap){
    //if(Verbosity()>0) 
    std::cout <<" xz cand out ev: " << _nevent << '[' << key << "] = " << value.first << " | " << value.second << std::endl;
  }

  for (const auto& [key, value] : yzouttrkmap){
    //if(Verbosity()>0) 
    std::cout <<" yz cand out ev: " << _nevent << '[' << key << "] = " << value.first << " | " << value.second << std::endl;
  }
  
  for (const auto& [xzkey, xzvalue] : xzouttrkmap){
    std::cout <<" xz cand out ev: " << _nevent << '[' << xzkey << "] = " << xzvalue.first << " | " << xzvalue.second << std::endl;

    //Loop over clusters and pick up clusters compatible with line parameters

    float xzint = -1*xzvalue.first;
    float xzsl = -1*xzvalue.second;
    
    for (const auto& [yzkey, yzvalue] : yzouttrkmap){
      std::cout <<" yyz cand out ev: " << _nevent << '[' << xzkey << "] = " << yzvalue.first << " | " << xzvalue.second << std::endl;
      
      //Loop over clusters and pick up clusters compatible with line parameters
      
      float yzint = -1*yzvalue.first;
      float yzsl = -1*yzvalue.second;
      
      // The shortest distance of a point from a circle is along the radial; line from the circle center to the point
      /*
      double xza = -xzsl;//-slope;
      double xzb = -1.0;
      double xzc = -tint;//-intercept;
      double xzr = 100;
      double xzx0 = -xza*xzc/(xza*xza+xzb*xzb), xzy0 = -xzb*xzc/(xza*xza+xzb*xzb);
      double xzax = NAN, xzay = NAN, xzbx = NAN, xzby = NAN;
      if ((xzc*xzc) > (xzr*xzr*(xza*xzxa+xzb*xzb)+0.00000001)){
	if(Verbosity()>0) cout << "no points " << endl;
      }
      else if (abs (xzc*xzc - xzr*xzr*(xza*xza+xzb*xzxb)) < 0.0000001) {
      if(Verbosity()>0){ puts ("1 point");
	cout << xzx0 << ' ' << xzy0 << endl;
      }
      }
      else {
	double xzd = xzr*xzr - xzc*xzc/(xza*xza+xzb*xzb);
	double xzmult = sqrt (xzd / (xza*xza+xzb*xzb));
	xzax = xzx0 + xzb * xzmult;
	xzbx = xzx0 - xzb * xzmult;
	xzay = xzy0 - xza * xzmult;
	xzby = xzy0 + xza * xzmult;
	if(Verbosity()>0){puts ("2 points");
	  cout << ax << ' ' << ay << '\n' << bx << ' ' << by << endl;
	}
      }

    float dist = 4;
    float alpha = atan(tsl);
    float offax = -1*dist*sin(alpha);
    float offay = -1*dist*cos(alpha);
    float offbx = 1*dist*sin(alpha);
    float offby = 1*dist*cos(alpha);
    vector<pointKey> trkclusters;
      */
      //    float trkdist = 2;
      /* float zstart = -110;
    float zend = 110;
    float xstart = xzint + zstart*xzsl; 
    float xend = xzint + zend*xzsl; 
    float ystart = yzint + zstart*yzsl; 
    float yend = yzint + zend*yzsl; */
    /*
    if(Verbosity()>0){ cout << ax << ' ' << ay << ' ' << bx << ' ' << by << endl;
      cout << "b1 = new TLine(" << ax +offax << ',' << ay +offay << ',' << bx -offbx << ',' << by - offby<<")" << endl;
      cout << "b2 = new TLine(" << ax -offax << ',' << ay -offay << ',' << bx +offbx << ',' << by + offby<<")" << endl;
      cout << "b1->Draw(\"same\")" << endl;
      cout << "b2->Draw(\"same\")" << endl;
      }*/
    vector<pointKey> trkclusters;
    vector<pointKey> lineclusters;
    rtree.query(bgi::intersects(box(point(-80,-80,-110),point(80,80,110))),std::back_inserter(lineclusters));
    // if(Verbosity()>0) 
    cout << "number line clus is " << lineclusters.size() << endl;
    
    for(vector<pointKey>::iterator clustertrk = lineclusters.begin();clustertrk!=lineclusters.end();++clustertrk){
      
      float ptx = clustertrk->first.get<0>();
      float pty = clustertrk->first.get<1>();
      float ptz = clustertrk->first.get<2>();
      
      float refx = xzint + ptz*xzsl;
      float refy = yzint + ptz*yzsl;
      float resx =  abs(ptx-refx);
      float resy =  abs(pty-refy);
      //      if(Verbosity()>0)
      if((resx<2)&&(resy<2)){
	cout << " x: " << ptx << " | " << refx << " y: " << pty << " | " << refy << " resx " << resx << " resy:" << resy << endl;
	trkclusters.push_back(*clustertrk);
      }
    }
    int nclus = 0;
    int nclus0 = 0;
    int nclus1 = 0;
    std::vector<float> dedxlist;

    cout << "number trk clus is " << trkclusters.size() << endl;
    if(trkclusters.size()>=_min_nclusters){
      cout << "setting keep true: " << trkclusters.size() << " > " << _min_nclusters << endl;
      for (const auto& trkclu : trkclusters){
	//      std::cout <<" yyz cand out ev: " << _nevent << '[' << xzkey << "] = " << yzvalue.first << " | " << xzvalue.second << std::endl;
	// for(vector<pointKey>::iterator clustertrk = trkclusters.begin();clustertrk!=trkclusters.end();++clustertrk){
	float tptx = trkclu.first.get<0>();
	float tpty = trkclu.first.get<1>();
	float tptz = trkclu.first.get<2>();
	//      float ptx = clustertrk->first.get<0>();
	//float pty = clustertrk->first.get<1>();
	// float ptz = clustertrk->first.get<2>();
	 nclus++;
	 if(tptz<0){
	   nclus0 = 0;
	 }
	 if(tptz>0){
	   nclus1 = 0;
	 }
	 float adc = trkclu.second;
	 
	 float thick = 1.098;
	 float r  = sqrt(tptx*tptx+tpty*tpty);
	 if(r>30&&r<=40.5){
	   thick = 0.57;
	 }else if(r>40.5&&r<60){
	   thick = 1.02;
	 } 
	 dedxlist.push_back(adc/thick);
	 //calc dedx
	 //_ntp_clutrk->Fill(_nevent,tptx , tpty, tptz);
      }
      sort(dedxlist.begin(), dedxlist.end());
      int trunc_min = 0;
      int trunc_max = (int)dedxlist.size()*0.7;
      float sumdedx = 0;
      int ndedx = 0;
      for(int j = trunc_min; j<=trunc_max;j++){
	sumdedx+=dedxlist.at(j);
	ndedx++;
      }
      sumdedx/=ndedx;
      TVector2 vec(xzint,yzint);
      // ev:nclus:nclus0:nclus1:r:phi:dedx:intxz:intyz:intrz:slxz:slyz:slrz
      float r = sqrt(xzint*xzint+ yzint*yzint);
      float xp = xzint+xzsl*100;
      float xn = xzint+xzsl*-100;
      float yp = yzint+yzsl*100;
      float yn = yzint+yzsl*-100;
      float rp = sqrt(xp*xp+yp*yp);
      float rn = sqrt(xn*xn+yn*yn);
      float phi = vec.Phi();
      float dr = rp-rn;
      
      _ntp_trk->Fill(_nevent,ntrack,nclus,nclus0,nclus1,r,phi,sumdedx,xzint,yzint,r,xzsl,yzsl,dr/200);

      for (const auto& trkclu : trkclusters){
	float tptx = trkclu.first.get<0>();
	float tpty = trkclu.first.get<1>();
	float tptz = trkclu.first.get<2>();
	float ftX[20] = {0};
	int nft = 0;
	ftX[nft++] = _nevent;
	ftX[nft++] =ntrack;
	ftX[nft++] =tptx;
	ftX[nft++] =tpty;
	ftX[nft++] = tptz;
	ftX[nft++] =nclus;
	ftX[nft++] =nclus0;
	ftX[nft++] =nclus1;
	ftX[nft++] =r;
	ftX[nft++] =phi;
	ftX[nft++] =sumdedx;
	ftX[nft++] =xzint;
	ftX[nft++] =yzint;
	ftX[nft++] =r;
	ftX[nft++] =xzsl;
	ftX[nft++] =yzsl;
	ftX[nft++] =dr/200;
	
	_ntp_clutrk->Fill(ftX);
      }
      ntrack++;
      //      keep_event = true;
    }
    /*
    //Assemble tracks
    if(_create_tracks){
      if(trkclusters.size()>=20){
	auto trackseed = std::make_unique<TrackSeed_v2>();
	for(const auto& cluskeys:trkclusters){
	  //      for(vector<pointKey>::iterator cluskeys =trkclusters.begin();cluskeys!=trkclusters.end();trkclusters++){
	  
	  trackseed->insert_cluster_key(cluskeys.second);
	}
	//    clean_chains.push_back(trackseed);
	// auto pseed = std::make_unique<TrackSeed_v2>(trackseed);
	m_seedContainer->insert(trackseed.get()); 
	cout << "number trk keys is " << trackseed->size_cluster_keys() << endl;
	numberofseeds++;
      }
    */
    }
  }

  
  cout << "number of seeds is " << ntrack << endl;

  /*  if(!keep_event){
    cout << " ABORT !keep " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  */
#endif
  return Fun4AllReturnCodes::EVENT_OK;
}


int PHStreakFinder::Setup(PHCompositeNode* topNode)
{
  cout << "Called Setup" << endl;
  cout << "topNode:" << topNode << endl;
  // PHTrackSeeding::Setup(topNode);
  //  int ret = GetNodes(topNode);
  //return ret;
  GetNodes(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHStreakFinder::End(PHCompositeNode* topNode)
{
  if(topNode==nullptr){
    std::cout << PHWHERE << "No topNode, bailing."
	      << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
    
    _tfile->cd();
    _ntp_clu->Write();
    _ntp_stub->Write();
    _ntp_trk->Write();
    _ntp_clutrk->Write();
    _tfile->Close();
    cout << "Called End " << endl;
 
  return Fun4AllReturnCodes::EVENT_OK;
}


