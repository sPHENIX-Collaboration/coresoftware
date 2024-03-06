/*!
 *  \file PHCosmicsFilter.cc
 
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
#include <trackbase_historic/TrackSeed_v1.h>

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

#include "PHCosmicsFilter.h"

#include "AssocInfoContainer.h"                         // for AssocInfoCont...

// trackbase_historic includes
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>  // for getLayer, clu...
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed_v1.h>


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
using pointKey = std::pair<point, TrkrDefs::cluskey>;
//typedef std::pair<point, TrkrDefs::cluskey> pointKey;

using myrtree = bgi::rtree<pointKey, bgi::quadratic<16>>;

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

vector<TrkrCluster*> clusterpoints;

PHCosmicsFilter::PHCosmicsFilter(const std::string& name)
  : SubsysReco(name)
{
  
}

int PHCosmicsFilter::GetNodes(PHCompositeNode* topNode)
{
  tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tGeometry)
  {
    std::cout << PHWHERE << "No acts reco geometry, bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    std::cout << PHWHERE << "No cluster container on node tree, bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHCosmicsFilter::Init(PHCompositeNode* topNode)
{
  if(topNode==nullptr){
    std::cout << PHWHERE << "No topNode, bailing."
	      << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _nevent = 0;
  if(_write_ntp){
  _tfile = new TFile("./costuple.root", "RECREATE");
  _ntp_cos = new TNtuple("ntp_cos", "cos event info","ev:x:y");
  _ntp_stub = new TNtuple("ntp_stub", "cos stub info","ev:int:sl:tan:x0");
  _ntp_max = new TNtuple("ntp_max", "cos stub info","ev:intmin:intmax:slmin:slmax");
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHCosmicsFilter::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }
  ret = createNodes(topNode);


  return ret;
}
int PHCosmicsFilter::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHCosmicsFilter::createNodes");
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

double PHCosmicsFilter::phiadd(double phi1, double phi2){
  double s = phi1+phi2;
  if(s>2*M_PI) {return s-2*M_PI;}
  else if(s<0) {return s+2*M_PI;}
  else {return s;}
}


double PHCosmicsFilter::phidiff(double phi1, double phi2){
  double d = phi1-phi2;
  if(d>M_PI) {return d-2*M_PI;}
  else if(d<-M_PI) {return d+2*M_PI;}
  else {return d;}
}

void PHCosmicsFilter::get_stub(const myrtree &search_rtree, float pointx, float pointy, int &count, double &slope, double &intercept){ //NOLINT
  float m1_dx = 4;
  float m1_dy = 4;
  vector<pointKey> boxclusters;
  search_rtree.query(bgi::intersects(box(point(pointx-m1_dx,pointy-m1_dy,-1),
				  point(pointx+m1_dx,pointy+m1_dy,1))),std::back_inserter(boxclusters));
  
  int nbox = boxclusters.size();
  
  int nhit = 0;
  double xsum = 0;
  double x2sum = 0;
  double ysum = 0;
  double xysum = 0;
  if(nbox>=2){
    for(vector<pointKey>::iterator pbox = boxclusters.begin();pbox!=boxclusters.end();++pbox){
      float boxx = pbox->first.get<0>();
      float boxy = pbox->first.get<1>();
      //      std::cout << "fprobe: " << count << " x:"  << boxx << " y: " << boxy << std::endl;
      nhit++;
      float r = boxx;
      float z = boxy;
      xsum = xsum + r;            // calculate sigma(xi)
      ysum = ysum + z;            // calculate sigma(yi)
      x2sum = x2sum + r*r;  // calculate sigma(x^2i)
      xysum = xysum + r * z;      // calculate sigma(xi*yi)
    }
  }
  
  const double denominator = ((x2sum * nhit) - (xsum*xsum));
  slope = (xysum * nhit - xsum * ysum) / denominator;   // calculate slope
  intercept = (x2sum * ysum - xsum * xysum) / denominator;  // calculate intercept
  count = nhit;
}

int PHCosmicsFilter::process_event(PHCompositeNode* topNode)
{
  if(topNode==nullptr){
    std::cout << PHWHERE << "No topNode, bailing."
	      << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
    
  _nevent++;
  std::vector<TrackSeed_v1> clean_chains;
  //Fill rtree
  bgi::rtree<pointKey, bgi::quadratic<16> > rtree;

  float slmin =  99999999999.9;
  float slmax = -99999999999.9;
  float intmin = 99999999999.9;
  float intmax =-99999999999.9;
  int num = 0;
  //  for(const auto& hitsetkey:_cluster_map->getHitSetKeys(TrkrDefs::TrkrId::tpcId)){
  if(_cluster_map->size()<_min_nclusters){
    std::cout << " not enough clusters in event: " << _cluster_map->size() << " < " << _min_nclusters << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  for(const auto& hitsetkey:_cluster_map->getHitSetKeys(TrkrDefs::TrkrId::tpcId)){
  //  for(const auto& hitsetkey:_cluster_map->getHitSetKeys()){
    auto range = _cluster_map->getClusters(hitsetkey);
    for( auto clusIter = range.first; clusIter != range.second; ++clusIter ){
      TrkrDefs::cluskey ckey = clusIter->first;
      TrkrCluster *cluster = clusIter->second;
      //      unsigned int layer = TrkrDefs::getLayer(ckey);
      /* if (layer < 7 || layer >= 55){
	if(Verbosity()>0) std::cout << "layer: " << layer << std::endl;
	continue;
	}*/
      // get global position, convert to Acts::Vector3 and store in map
      const Acts::Vector3 globalpos_d =  tGeometry->getGlobalPosition(ckey, cluster);
      
      const Acts::Vector3 globalpos = { globalpos_d.x(), globalpos_d.y(), globalpos_d.z()};
      //      if(layer>=51&&layer<=55)
      num++;
      //      std::cout << "#: " << num << " l xyz " << layer << " | " << globalpos_d.x() << " | " << globalpos_d.y() << " | " <<  globalpos_d.z() <<std::endl;
      if(_write_ntp)_ntp_cos->Fill(_nevent,globalpos_d.x(), globalpos_d.y());
      /*
      const double clus_phi = get_phi( globalpos );      
      const double clus_eta = get_eta( globalpos );
      const double clus_l = layer;  
      */

      std::vector<pointKey> testduplicate;
      rtree.query(bgi::intersects(box(point(globalpos_d.x()-0.001,globalpos_d.y()-0.001,-1),
				      point(globalpos_d.x()+0.001,globalpos_d.y()+0.001,1))),std::back_inserter(testduplicate));
      if (!testduplicate.empty()){
	continue;
      }

      rtree.insert(std::make_pair(point(globalpos.x() , globalpos.y(), 0.0), ckey));
      
    }
  }

  //Get all clusters from rtree, fit stubs around clusters, fill stub tree
  vector<pointKey> allclusters;
  
  rtree.query(bgi::intersects(box(point(-80,-80,-1),point(80,80,1))),std::back_inserter(allclusters));
  if(Verbosity()>0) cout << "number clus is " << allclusters.size() << endl;

  bgi::rtree<pointKey, bgi::quadratic<16> > rtree_stub;

  for(vector<pointKey>::iterator cluster = allclusters.begin();cluster!=allclusters.end();++cluster){
    float pointx = cluster->first.get<0>();
    float pointy = cluster->first.get<1>();
    int fcount = 0;
    double fslope = 0;
    double fintercept = 0;
    //calc slope and intersect from 5 cluster stubs
    get_stub(rtree, pointx, pointy, fcount, fslope, fintercept);
    if(finite(fslope)){
	rtree_stub.insert(std::make_pair(point(fintercept , fslope, 0.0), 0));
	float tana = atan(fslope);
	if(_write_ntp)_ntp_stub->Fill(_nevent,fintercept,fslope,tana,fintercept);
	if(fslope>slmax)slmax = fslope;
	if(fslope<slmin)slmin = fslope;
	if(fintercept>intmax)intmax = fintercept;
	if(fintercept<intmin)intmin = fintercept;
      }
  }


  
  //  int out = 0;
  // find clusters of slope intectept pairs
  map <int,pair<float,float>> outtrkmap;
  int nout = 0;

  while(rtree_stub.size()>10){
    std::vector<pointKey> allstubs;
    rtree_stub.query(bgi::intersects(box(point(intmin,slmin,-1),point(intmax,slmax,1))),std::back_inserter(allstubs));
    map <int,pair<float,float>> trkmap;

    float int_width = (intmax - intmin)/20;
    float sl_width = (slmax - slmin)/20;
    for(vector<pointKey>::iterator stub = allstubs.begin();stub!=allstubs.end();++stub){
      float pint = stub->first.get<0>();
      float psl = stub->first.get<1>();
      vector<pointKey> trkcand;
      
      rtree_stub.query(bgi::intersects(box(point(pint-int_width,psl-sl_width,-1),point(pint+int_width,psl+sl_width,1))),std::back_inserter(trkcand));
      
      int ntrk = trkcand.size();
      int count = 0;
      float intsum = 0;
      float slsum = 0;
      if(ntrk>=5){
	for(vector<pointKey>::iterator ptrk = trkcand.begin();ptrk!=trkcand.end();++ptrk){
	  float trkint = ptrk->first.get<0>();
	  float trksl = ptrk->first.get<1>();
	  if(Verbosity()>0) cout<< "    stub " << ntrk
				<< " int: " << trkint 
				<< " sl: "  << trksl 
				<< endl;
	  
	  intsum = intsum + trkint;            // calculate sigma(xi)
	  slsum = slsum + trksl;            // calculate sigma(yi)
	  count++;
	}
	float mint = (intsum/count);
	float msl  = (slsum/count);
	trkmap[ntrk] = std::make_pair(mint,msl); 
      
	if(Verbosity()>0) cout<< " stub in box " << ntrk
			      << " int: " << pint 
			      << " sl: " << psl 
			      << " intwi: " << int_width 
			      << " slwi: " << sl_width
			      << " mint " << mint
			      << " msl " << msl
			      << endl;
      }
    }
    
    if( trkmap.size()>0){
      int size_before = rtree_stub.size();
      if(Verbosity()>0) cout << "mapend: " << trkmap.rbegin()->first << " | " << (trkmap.rbegin()->second).first << " | "  << trkmap.rbegin()->second.second << endl;
      // for (const auto& [key, value] : trkmap)
      //	std::cout <<" ev: " << _nevent << '[' << key << "] = " << value.first << " | " << value.second << std::endl;
      //remove stubs in box
      float rmint = trkmap.rbegin()->second.first;
      float rmsl = trkmap.rbegin()->second.second;
      vector<pointKey> rmcand;
      rtree_stub.query(bgi::intersects(box(point(rmint-int_width,rmsl-sl_width,-1),point(rmint+int_width,rmsl+sl_width,1))),std::back_inserter(rmcand));
      for(vector<pointKey>::iterator rmstub = rmcand.begin();rmstub!=rmcand.end();++rmstub){
	float rmpint = rmstub->first.get<0>();
	float rmpsl = rmstub->first.get<1>();
	if(Verbosity()>0) cout<< "    rm " <<  " int: " << rmpint 
			      << " sl: " << rmpsl 
			      << endl;

	//auto rmpoint = std::make_pair(point(rmpint , rmpsl, 0.0), 0);
	//rtree_stub.remove(rmpoint);
	//	pointKey match;
	rtree_stub.remove(*rmstub);
	/*	if (!rtree_stub.query(bgi::nearest(rmpoint, 1), &match))
	  continue;
	*/
	/*
	if (bg::distance(rmpoint, match.first) <= 0.5) {
	  rtree_stub.remove(match);
	 
	  */
      }
      int size_after = rtree_stub.size();
      if(size_before == size_after) {break;}
      outtrkmap[nout++] = std::make_pair(rmint,rmsl);
      cout << " tree sixe after remove: " << rtree_stub.size() << endl;
    }else{
      break;
    }
    trkmap.clear();
  }

  int numberofseeds = 0;
  bool keep_event = false;
  for (const auto& [key, value] : outtrkmap){
    if(Verbosity()>0) std::cout <<" out ev: " << _nevent << '[' << key << "] = " << value.first << " | " << value.second << std::endl;
    //Loop over clusters and pick up clusters compatible with line parameters

    float tint = -1*value.first;
    float tsl = -1*value.second;
    
    //  float liney = tsl*ptx + tint;
    
    // double ry = pty;
    // double zx = ptx;
    
    // The shortest distance of a point from a circle is along the radial; line from the circle center to the point
    
    double a = -tsl;//-slope;
    double b = -1.0;
    double c = -tint;//-intercept;
    double r = 100;
    double x0 = -a*c/(a*a+b*b), y0 = -b*c/(a*a+b*b);
    double ax = NAN, ay = NAN, bx = NAN, by = NAN;
    if ((c*c) > (r*r*(a*a+b*b)+0.00000001)){
      if(Verbosity()>0) cout << "no points " << endl;
    }
    else if (abs (c*c - r*r*(a*a+b*b)) < 0.0000001) {
      if(Verbosity()>0){ puts ("1 point");
	cout << x0 << ' ' << y0 << endl;
      }
    }
    else {
      double d = r*r - c*c/(a*a+b*b);
      double mult = sqrt (d / (a*a+b*b));
      ax = x0 + b * mult;
      bx = x0 - b * mult;
      ay = y0 - a * mult;
      by = y0 + a * mult;
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
    
    //float trkdist = 4;
    vector<pointKey> lineclusters;
    if(Verbosity()>0){ cout << ax << ' ' << ay << ' ' << bx << ' ' << by << endl;
      cout << "b1 = new TLine(" << ax +offax << ',' << ay +offay << ',' << bx -offbx << ',' << by - offby<<")" << endl;
      cout << "b2 = new TLine(" << ax -offax << ',' << ay -offay << ',' << bx +offbx << ',' << by + offby<<")" << endl;
      cout << "b1->Draw(\"same\")" << endl;
      cout << "b2->Draw(\"same\")" << endl;
    }
    rtree.query(bgi::intersects(box(point(-80,-80,-1),point(80,80,1))),std::back_inserter(lineclusters));
    if(Verbosity()>0)  cout << "number line clus is " << lineclusters.size() << endl;
    
    // Check if track hits MVTX (dist to 0,0 < 2)

    float dist_origin = std::abs(c)/sqrt(a*a+b*b);
    if(_max_dist_to_origin>0){
      if(dist_origin>_max_dist_to_origin){
	cout << "dist continue: " << dist_origin << " > " << _max_dist_to_origin << endl;
	continue;
      }
    }
    for(vector<pointKey>::iterator clustertrk = lineclusters.begin();clustertrk!=lineclusters.end();++clustertrk){
      
      float ptx = clustertrk->first.get<0>();
      float pty = clustertrk->first.get<1>();
      

      float res =  std::abs(a*ptx+b*pty+c)/sqrt(a*a+b*b);
      //float res = std::abs(tsl*ptx+(-1*pty)+tint)/sqrt((tsl*tsl+1));
      if(Verbosity()>0) cout << " x: " << ptx << " y: " << pty << " res " << res << endl;
      if(res<4){
	trkclusters.push_back(*clustertrk);
      }
    }
    cout << "number trk clus is " << trkclusters.size() << endl;
    if(trkclusters.size()>=_min_nclusters){
      cout << "setting keep true: " << trkclusters.size() << " > " << _min_nclusters << endl;
      keep_event = true;
    }

    //Assemble tracks
    if(_create_tracks){
      if(trkclusters.size()>=20){
	auto trackseed = std::make_unique<TrackSeed_v1>();
	for(const auto& cluskeys:trkclusters){
	  //      for(vector<pointKey>::iterator cluskeys =trkclusters.begin();cluskeys!=trkclusters.end();trkclusters++){
	  
	  trackseed->insert_cluster_key(cluskeys.second);
	}
	//    clean_chains.push_back(trackseed);
	// auto pseed = std::make_unique<TrackSeed_v1>(trackseed);
	m_seedContainer->insert(trackseed.get()); 
	cout << "number trk keys is " << trackseed->size_cluster_keys() << endl;
	numberofseeds++;
      }
    }
  }
  
  cout << "number of seeds is " << numberofseeds << endl;
  if(_write_ntp){
    _ntp_max->Fill(_nevent,intmin,intmax,slmin,slmax);
  }
  if(!keep_event){
    cout << " ABORT !keep " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}


int PHCosmicsFilter::Setup(PHCompositeNode* topNode)
{
  cout << "Called Setup" << endl;
  cout << "topNode:" << topNode << endl;
  // PHTrackSeeding::Setup(topNode);
  //  int ret = GetNodes(topNode);
  //return ret;
  GetNodes(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCosmicsFilter::End(PHCompositeNode* topNode)
{
  if(topNode==nullptr){
    std::cout << PHWHERE << "No topNode, bailing."
	      << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
    
  if(_write_ntp){
    _tfile->cd();
    _ntp_cos->Write();
    _ntp_stub->Write();
    _ntp_max->Write();
    _tfile->Close();
    cout << "Called End " << endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}


