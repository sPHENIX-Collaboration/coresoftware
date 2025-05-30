/*!
 *  \file PHLineLaserReco.cc
 
 *  \author Christof Roland 
 */


//begin

#include "PHLineLaserReco.h"

#include "AssocInfoContainer.h"                         // for AssocInfoCont...

// sPHENIX Geant4 includes
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/TrkrDefs.h>  // for getLayer, clu...
#include <trackbase/LaserCluster.h>
#include <trackbase/LaserClusterContainer.h>
#include <trackbase/LaserClusterContainerv1.h>
#include <trackbase/LaserClusterv1.h>

#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>

// sPHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHTimer.h>                              // for PHTimer
#include <phool/getClass.h>

//ROOT includes for debugging
#include <TFile.h>
#include <TH1.h>
#include <TNtuple.h>
#include <TAxis.h>
#include <TGraph.h>

#include <Eigen/Core>                  // for Matrix
#include <Eigen/Dense>

//BOOST for combi seeding
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

// standard includes
#include <algorithm>
#include <cfloat>
#include <climits>                                     // for UINT_MAX
#include <cmath>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <memory>
#include <set>                                          // for set
#include <string>                                  // for string
#include <tuple>
#include <utility>
#include <vector>
                                   // for pair, make_pair

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


using point = bg::model::point<float, 3, bg::cs::cartesian>;
//typedef bg::model::point<float, 3, bg::cs::cartesian> point;
using box = bg::model::box<point>;
//typedef bg::model::box<point> box;
using pointInd = std::pair<point, TrkrDefs::cluskey>;
//typedef std::pair<point, TrkrDefs::cluskey> pointKey;

using myrtree = bgi::rtree<pointInd, bgi::quadratic<16>>;


#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp


using namespace std;

vector<LaserCluster*> laserpoints;

PHLineLaserReco::PHLineLaserReco(const std::string& name)
  : SubsysReco(name)
{
  
}

int PHLineLaserReco::GetNodes(PHCompositeNode* topNode)
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
int PHLineLaserReco::Init(PHCompositeNode* topNode)
{
  if(topNode==nullptr){
    std::cout << PHWHERE << "No topNode, bailing."
	      << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _nevent = 0;
  if(_write_ntp){
  _tfile = new TFile("./costuple.root", "RECREATE");
  _ntp_cos = new TNtuple("ntp_cos", "cos event info","ev:x:y:z:adc:maxadc:size:lsize:phisize:tsize");
  _ntp_stub = new TNtuple("ntp_stub", "cos stub info","ev:int:sl:tan:x0:count");
  _ntp_max = new TNtuple("ntp_max", "cos stub info","ev:intmin:intmax:slmin:slmax");
  _ntp_trk = new TNtuple("ntp_trk", "laser track info","ev:int:sl:intfit:slfit:nclus");
  //  _ntp_trk_clus = new TNtuple("ntp_trk_clus", "laser clus","ev:x:y:z:adc:maxadc:size:lsize:phisize:tsize:intfit:slfit:nclus:zfirst:zlast");
  _ntp_trk_hit = new TNtuple("ntp_trk_hit", "laser hitinfo","ev:x:y:z:adc:lay:phibin:zbin:intfit:slfit:nclus:zfirst:zlast");
  }

  m_hittree = new TTree("hittree", "A tree with all hits");
  m_hittree->Branch("event", &m_nevent, "m_nevent/I");
  m_hittree->Branch("x", &m_hitx, "m_hitx/F");
  m_hittree->Branch("y", &m_hity, "m_hity/F");
  m_hittree->Branch("z", &m_hitz, "m_hitz/F");
  m_hittree->Branch("layer", &m_hitlayer, "m_hitlayer/I");
  m_hittree->Branch("adc", &m_hitadc, "m_hitadc/I");
  m_hittree->Branch("pad", &m_hitpad, "m_hitpad/I");
  m_hittree->Branch("tbin", &m_hittbin, "m_hittbin/I");
  m_hittree->Branch("slopexy", &m_slopexy, "m_slopexy/F");
  m_hittree->Branch("interxy", &m_interxy, "m_interxy/F");
  m_hittree->Branch("slopexz", &m_slopexz, "m_slopexz/F");
  m_hittree->Branch("interxz", &m_interxz, "m_interxz/F");
  m_hittree->Branch("slopeyz", &m_slopeyz, "m_slopeyz/F");
  m_hittree->Branch("interyz", &m_interyz, "m_interyz/F");
  m_hittree->Branch("nclus", &m_nclus, "m_nclus/I");
  m_hittree->Branch("zlast", &m_zlast, "m_zlast/F");
  m_hittree->Branch("zfirst", &m_zfirst, "m_first/F");

  m_clustree = new TTree("clustree", "A tree with all hits");
  m_clustree->Branch("event", &m_nevent, "m_nevent/I");
  m_clustree->Branch("x", &m_clux, "m_clux/F");
  m_clustree->Branch("y", &m_cluy, "m_cluy/F");
  m_clustree->Branch("z", &m_cluz, "m_cluz/F");
  m_clustree->Branch("adc", &m_cluadc, "m_cluadc/I");
  m_clustree->Branch("maxadc", &m_clumaxadc, "m_clumaxadc/I");
  m_clustree->Branch("size", &m_size, "m_size/I");
  m_clustree->Branch("lsize", &m_sizel, "m_sizel/I");
  m_clustree->Branch("phisize", &m_sizephi, "m_sizephi/I");
  m_clustree->Branch("tsize", &m_sizet, "m_sizet/I");
  m_clustree->Branch("slopexy", &m_slopexy, "m_slopexy/F");
  m_clustree->Branch("interxy", &m_interxy, "m_interxy/F");
  m_clustree->Branch("slopexz", &m_slopexz, "m_slopexz/F");
  m_clustree->Branch("interxz", &m_interxz, "m_interxz/F");
  m_clustree->Branch("slopeyz", &m_slopeyz, "m_slopeyz/F");
  m_clustree->Branch("interyz", &m_interyz, "m_interyz/F");
  m_clustree->Branch("nclus", &m_nclus, "m_nclus/I");
  m_clustree->Branch("zlast", &m_zlast, "m_zlast/F");
  m_clustree->Branch("zfirst", &m_zfirst, "m_first/F");
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHLineLaserReco::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }
  ret = createNodes(topNode);


  return ret;
}
int PHLineLaserReco::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHLineLaserReco::createNodes");
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

double PHLineLaserReco::phiadd(double phi1, double phi2){
  double s = phi1+phi2;
  if(s>2*M_PI) {return s-2*M_PI;}
  else if(s<0) {return s+2*M_PI;}
  else {return s;}
}


double PHLineLaserReco::phidiff(double phi1, double phi2){
  double d = phi1-phi2;
  if(d>M_PI) {return d-2*M_PI;}
  else if(d<-M_PI) {return d+2*M_PI;}
  else {return d;}
}

void PHLineLaserReco::get_stub(const myrtree &search_rtree, float pointx, float pointy, int &count, double &slope, double &intercept){ //NOLINT
  float m1_dx = 6;
  float m1_dy = 6;
  vector<pointInd> boxclusters;
  search_rtree.query(bgi::intersects(box(point(pointx-m1_dx,pointy-m1_dy,-1),
				  point(pointx+m1_dx,pointy+m1_dy,1))),std::back_inserter(boxclusters));
  
  int nbox = boxclusters.size();
  
  int nhit = 0;
  double xsum = 0;
  double x2sum = 0;
  double ysum = 0;
  double xysum = 0;
  if(nbox>=2){
    for(vector<pointInd>::iterator pbox = boxclusters.begin();pbox!=boxclusters.end();++pbox){
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

int PHLineLaserReco::process_event(PHCompositeNode* topNode)
{
  std::cout << "start processing event: " << std::endl;
  if(topNode==nullptr){
    std::cout << PHWHERE << "No topNode, bailing."
	      << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
    
  _nevent++;
  std::vector<TrackSeed_v2> clean_chains;
  //Fill rtree
  bgi::rtree<pointInd, bgi::quadratic<16> > rtree;

  float slmin =  99999999999.9;
  float slmax = -99999999999.9;
  float intmin = 99999999999.9;
  float intmax =-99999999999.9;
//  int num = 0;
  //  for(const auto& hitsetkey:_cluster_map->getHitSetKeys(TrkrDefs::TrkrId::tpcId)){
  if(_cluster_map->size()<_min_nclusters){
    std::cout << " not enough clusters in event: " << _cluster_map->size() << " < " << _min_nclusters << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto range = _cluster_map->getClusters();
  std::cout << "_nevent " << _nevent << "n clusters" << std::distance(range.first,range.second) << std::endl; 
  for( auto clusIter = range.first; clusIter != range.second; ++clusIter ){
    LaserCluster *cluster = clusIter->second;
    TrkrDefs::cluskey lckey = clusIter->first;
    const Acts::Vector3 globalpos_d =  {cluster->getX() , cluster->getY(), cluster->getZ()};
    const Acts::Vector3 globalpos = { globalpos_d.x(), globalpos_d.y(), globalpos_d.z()};
    //    ev:x:y:z:adc:maxadc:size:lsize:phisize:tsize
    //    if(_write_ntp)_ntp_cos->Fill(_nevent,globalpos_d.x(), globalpos_d.y(),globalpos_d.z(),);
    /*
      const double clus_phi = get_phi( globalpos );      
      const double clus_eta = get_eta( globalpos );
      const double clus_l = layer;  
    */
    //    int side = -1;
    unsigned int adcSum = 0.0;
    unsigned int maxAdc = 0.0;
    int iphimin = 6666, iphimax = -1;
    int ilaymin = 6666, ilaymax = -1;
    float itmin = 66666666.6, itmax = -6666666666.6;
    unsigned int nHits = cluster->getNhits();

    for(unsigned int i  = 0; i< nHits;i++){
      //      float coords[3] = {cluster->getHitX(i), cluster->getHitY(i), cluster->getHitZ(i)};
      int iphi = cluster->getHitIPhi(i);
      int it = cluster->getHitIT(i);
      int ilay = cluster->getHitLayer(i);
      unsigned int adc =  cluster->getHitAdc(i);
      if(iphi<iphimin){iphimin = iphi;}
      if(iphi>iphimax){iphimax = iphi;}
      if(ilay<ilaymin){ilaymin = ilay;}
      if(ilay>ilaymax){ilaymax = ilay;}
      if(it<itmin){itmin = it;}
      if(it>itmax){itmax = it;}
      adcSum += adc;
      if(adc>maxAdc){
	maxAdc = adc;
      }
    }
    int phisize =  iphimax - iphimin + 1;
    int lsize =  ilaymax - ilaymin + 1;
    int tsize = itmax - itmin +1;
    
    if(_write_ntp)
      //    ev:x:y:z:adc:maxadc:size:lsize:phisize:tsize
      _ntp_cos->Fill(_nevent,globalpos_d.x(), globalpos_d.y(),globalpos_d.z(),adcSum,maxAdc,nHits,lsize,phisize,tsize);
    if(lsize>2||maxAdc<50||tsize<=1)
      continue;
    std::vector<pointInd> testduplicate;
    rtree.query(bgi::intersects(box(point(globalpos_d.x()-0.001,globalpos_d.y()-0.001,-1),
				    point(globalpos_d.x()+0.001,globalpos_d.y()+0.001,1))),std::back_inserter(testduplicate));
    if (!testduplicate.empty()){
      continue;
    }
    rtree.insert(std::make_pair(point(globalpos.x() , globalpos.y(), 0.0), lckey)); 
  }

  //Get all clusters from rtree, fit stubs around clusters, fill stub tree
  vector<pointInd> allclusters;
  
  rtree.query(bgi::intersects(box(point(-80,-80,-1),point(80,80,1))),std::back_inserter(allclusters));
  if(Verbosity()>0) cout << "number clus is " << allclusters.size() << endl;

  bgi::rtree<pointInd, bgi::quadratic<16> > rtree_stub;

  for(vector<pointInd>::iterator cluster = allclusters.begin();cluster!=allclusters.end();++cluster){
    float pointx = cluster->first.get<0>();
    float pointy = cluster->first.get<1>();
    int fcount = 0;
    double fslope = 0;
    double fintercept = 0;
    //calc slope and intersect from 5 cluster stubs
    get_stub(rtree, pointx, pointy, fcount, fslope, fintercept);
    if(finite(fslope)&&fcount>5){
	rtree_stub.insert(std::make_pair(point(fintercept , fslope, 0.0), 0));
	float tana = atan(fslope);
	if(_write_ntp)
	  _ntp_stub->Fill(_nevent,fintercept,fslope,tana,fintercept,fcount);
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

  while(rtree_stub.size()>3){
    std::vector<pointInd> allstubs;
    rtree_stub.query(bgi::intersects(box(point(intmin,slmin,-1),point(intmax,slmax,1))),std::back_inserter(allstubs));
    map <int,pair<float,float>> trkmap;

    float int_width = (intmax - intmin)/10;
    float sl_width = (slmax - slmin)/10;
    for(vector<pointInd>::iterator stub = allstubs.begin();stub!=allstubs.end();++stub){
      float pint = stub->first.get<0>();
      float psl = stub->first.get<1>();
      vector<pointInd> trkcand;
      
      rtree_stub.query(bgi::intersects(box(point(pint-int_width,psl-sl_width,-1),point(pint+int_width,psl+sl_width,1))),std::back_inserter(trkcand));
      
      int ntrk = trkcand.size();
      int count = 0;
      float intsum = 0;
      float slsum = 0;
      if(ntrk>=5){
	for(vector<pointInd>::iterator ptrk = trkcand.begin();ptrk!=trkcand.end();++ptrk){
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
      vector<pointInd> rmcand;
      rtree_stub.query(bgi::intersects(box(point(rmint-int_width,rmsl-sl_width,-1),point(rmint+int_width,rmsl+sl_width,1))),std::back_inserter(rmcand));
      for(vector<pointInd>::iterator rmstub = rmcand.begin();rmstub!=rmcand.end();++rmstub){
	float rmpint = rmstub->first.get<0>();
	float rmpsl = rmstub->first.get<1>();
	if(Verbosity()>0) cout<< "    rm " <<  " int: " << rmpint 
			      << " sl: " << rmpsl 
			      << endl;

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
      outtrkmap.insert(std::make_pair(nout,std::make_pair(rmint,rmsl)));
      nout++;
      cout << " tree sixe after remove: " << rtree_stub.size() << endl;
    }else{
      break;
    }
    trkmap.clear();
  }

  int numberofseeds = 0;
  //  bool keep_event = false;
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
    
    //float trkdist = 4;
    vector<pointInd> lineclusters;
    if(Verbosity()>0){ cout << ax << ' ' << ay << ' ' << bx << ' ' << by << endl;
      cout << "auto b0 = new TLine(" << ax << ',' << ay << ',' << bx  << ',' << by <<")" << endl;
      cout << "auto b1 = new TLine(" << ax +offax << ',' << ay +offay << ',' << bx -offbx << ',' << by - offby<<")" << endl;
      cout << "auto b2 = new TLine(" << ax -offax << ',' << ay -offay << ',' << bx +offbx << ',' << by + offby<<")" << endl;
      cout << "b0->Draw(\"same\")" << endl;
      cout << "b1->Draw(\"same\")" << endl;
      cout << "b2->Draw(\"same\")" << endl;
    }
    rtree.query(bgi::intersects(box(point(-80,-80,-1),point(80,80,1))),std::back_inserter(lineclusters));
    if(Verbosity()>0)  cout << "number line clus is " << lineclusters.size() << endl;
    
    // Check if track hits MVTX (dist to 0,0 < 2)

    //first fit
    int nhitff = 0;
    double xsumff = 0;
    double ysumff = 0;
    double zsumff = 0;
    double x2sumff = 0;
    double y2sumff = 0;
    double z2sumff = 0;
    double xysumff = 0;
    double xzsumff = 0;
    double yzsumff = 0;
    double zfirst = 66666666;
    double zlast = -66666666;
    vector<pointInd> trkclusters;
    unsigned int trk_nclus = 0;

    //    fitstraight lines first
    for(vector<pointInd>::iterator clustertrk = lineclusters.begin();clustertrk!=lineclusters.end();++clustertrk){
      
      float ptx = clustertrk->first.get<0>();
      float pty = clustertrk->first.get<1>();
      TrkrDefs::cluskey tcind = clustertrk->second;
      LaserCluster *las_clus = _cluster_map->findCluster(tcind);      
      float ptz = las_clus->getZ();
      
      float res =  std::abs(a*ptx+b*pty+c)/sqrt(a*a+b*b);
      //float res = std::abs(tsl*ptx+(-1*pty)+tint)/sqrt((tsl*tsl+1));
      // if(Verbosity()>0)
      cout << " x: " << ptx << " y: " << pty << " res " << res << endl;
      if(res<4){
	float boxx = ptx;
	float boxy = pty;
	float boxz = ptz;
        std::cout << "fprobe x:"  << boxx << " y: " << boxy <<  " z: " << boxz << std::endl;
	nhitff++;
	float xff = boxx;
	float yff = boxy;
	float zff = boxz;

	xsumff = xsumff + xff;            // calculate sigma(xi)
	ysumff = ysumff + yff;            // calculate sigma(yi)
	zsumff = zsumff + zff;            // calculate sigma(yi)

	x2sumff = x2sumff + xff*xff;      // calculate sigma(x^2i)
	y2sumff = y2sumff + yff*yff;      // calculate sigma(x^2i)
       	z2sumff = z2sumff + zff*zff;      // calculate sigma(x^2i)

	xysumff = xysumff + xff * yff;    // calculate sigma(xi*yi)
	xzsumff = xzsumff + xff * zff;    // calculate sigma(xi*yi)
	yzsumff = yzsumff + yff * zff;    // calculate sigma(xi*yi)

	//	trkclusters.push_back(*clustertrk);
	//trk_nclus++;
	if(ptz<zfirst){zfirst = ptz;}
	if(ptz<zlast){zlast = ptz;}
      }
    }
    cout << "rescheck done " <<  endl;
    cout << "number trk clus is " << trkclusters.size() << " " << trk_nclus<< endl;
  


    double denominatorxff = ((x2sumff * nhitff) - (xsumff*xsumff));
    double slopexyff = (xysumff* nhitff - xsumff * ysumff) / denominatorxff;   // calculate slope
    double interceptxyff = (x2sumff * ysumff - xsumff * xysumff) / denominatorxff;  // calculate intercept

    double denominatorxzff = ((z2sumff * nhitff) - (zsumff*zsumff));
    double slopexzff = (xzsumff* nhitff - zsumff * xsumff) / denominatorxzff;   // calculate slope
    double interceptxzff = (z2sumff * xsumff - zsumff * xzsumff) / denominatorxzff;  // calculate intercept

    double denominatoryzff = ((z2sumff * nhitff) - (zsumff*zsumff));
    double slopeyzff = (yzsumff* nhitff - zsumff * ysumff) / denominatoryzff;   // calculate slope
    double interceptyzff = (z2sumff * ysumff - zsumff * yzsumff) / denominatoryzff;  // calculate intercept

    _ntp_trk->Fill(_nevent,tint,tsl,interceptxyff,slopexyff,trk_nclus);

    //  _ntp_trk_clus = new TNtuple("ntp_trk_clus", "laser clus","ev:x:y:z:adc:maxadc:size:lsize:phisize:tsize:intfit:slfit:nclus:zfirst:zlast");
    //_ntp_trk_hit = new TNtuple("ntp_trk_hit", "laser hitinfo","ev:x:y:z:adc:maxadc:lay:phibin:zbin:intfit:slfit:nclus:zfirst:zlast");

    //clean based on z position
    xsumff = 0;
    ysumff = 0;
    zsumff = 0;
    x2sumff = 0;
    y2sumff = 0;
    z2sumff = 0;
    xysumff = 0;
    xzsumff = 0;
    yzsumff = 0;

    // The shortest distance of a point from a circle is along the radial; line from the circle center to the point
    
    double axz =  slopexzff;//-slope;
    double bxz = -1.0;
    double cxz =  interceptxzff;//-intercept;

    double ayz =  slopeyzff;//-slope;
    double byz = -1.0;
    double cyz =  interceptyzff;//-intercept;

    for(vector<pointInd>::iterator clustertrk = lineclusters.begin();clustertrk!=lineclusters.end();++clustertrk){
      
      float ptx = clustertrk->first.get<0>();
      float pty = clustertrk->first.get<1>();
      TrkrDefs::cluskey tcind = clustertrk->second;
      LaserCluster *las_clus = _cluster_map->findCluster(tcind);      
      float ptz = las_clus->getZ();
      
      float res_xz =  std::abs(axz*ptz+bxz*ptx+cxz)/sqrt(axz*axz+bxz*bxz);
      float res_yz =  std::abs(ayz*ptz+byz*pty+cyz)/sqrt(ayz*ayz+byz*byz);
      //float res = std::abs(tsl*ptx+(-1*pty)+tint)/sqrt((tsl*tsl+1));
      // if(Verbosity()>0)
      cout << " x: " << ptx << " z: " << ptz << " res_xz " << res_xz << endl;
      cout << " y: " << pty << " z: " << ptz << " res_yz " << res_yz << endl;
      if(res_yz<4&&res_xz<4){
	float boxx = ptx;
	float boxy = pty;
	float boxz = ptz;
        std::cout << "fprobe x:"  << boxx << " y: " << boxy <<  " z: " << boxz << std::endl;
	nhitff++;
	float xff = boxx;
	float yff = boxy;
	float zff = boxz;

	xsumff = xsumff + xff;            // calculate sigma(xi)
	ysumff = ysumff + yff;            // calculate sigma(yi)
	zsumff = zsumff + zff;            // calculate sigma(yi)

	x2sumff = x2sumff + xff*xff;      // calculate sigma(x^2i)
	y2sumff = y2sumff + yff*yff;      // calculate sigma(x^2i)
       	z2sumff = z2sumff + zff*zff;      // calculate sigma(x^2i)

	xysumff = xysumff + xff * yff;    // calculate sigma(xi*yi)
	xzsumff = xzsumff + xff * zff;    // calculate sigma(xi*yi)
	yzsumff = yzsumff + yff * zff;    // calculate sigma(xi*yi)

	trkclusters.push_back(*clustertrk);
	trk_nclus++;
	if(ptz<zfirst){zfirst = ptz;}
	if(ptz<zlast){zlast = ptz;}
      }
    }
    cout << "rescheck done " <<  endl;
    cout << "number trk clus is " << trkclusters.size() << " " << trk_nclus<< endl;
  
    denominatorxff = ((x2sumff * nhitff) - (xsumff*xsumff));
    slopexyff = (xysumff* nhitff - xsumff * ysumff) / denominatorxff;   // calculate slope
    interceptxyff = (x2sumff * ysumff - xsumff * xysumff) / denominatorxff;  // calculate intercept

    denominatorxzff = ((z2sumff * nhitff) - (zsumff*zsumff));
    slopexzff = (xzsumff* nhitff - xsumff * zsumff) / denominatorxzff;   // calculate slope
    interceptxzff = (z2sumff * xsumff - zsumff * xzsumff) / denominatorxzff;  // calculate intercept

    denominatoryzff = ((z2sumff * nhitff) - (zsumff*zsumff));
    slopeyzff = (yzsumff* nhitff - zsumff * ysumff) / denominatoryzff;   // calculate slope
    interceptyzff = (z2sumff * ysumff - zsumff * yzsumff) / denominatoryzff;  // calculate intercept

    _ntp_trk->Fill(_nevent,tint,tsl,interceptxyff,slopexyff,trk_nclus);

    //  _ntp_trk_clus = new TNtuple("ntp_trk_clus", "laser clus","ev:x:y:z:adc:maxadc:size:lsize:phisize:tsize:intfit:slfit:nclus:zfirst:zlast");
    //_ntp_trk_hit = new TNtuple("ntp_trk_hit", "laser hitinfo","ev:x:y:z:adc:maxadc:lay:phibin:zbin:intfit:slfit:nclus:zfirst:zlast"
    
    for(vector<pointInd>::iterator trkclusiter = trkclusters.begin();trkclusiter!=trkclusters.end();++trkclusiter){
      float tcx = trkclusiter->first.get<0>();
      float tcy = trkclusiter->first.get<1>();
      float tcz = trkclusiter->first.get<2>();
      TrkrDefs::cluskey tcind = trkclusiter->second;
      LaserCluster *las_clus = _cluster_map->findCluster(tcind);      
      //int nth = 0;
      //int ntc = 0;
      // float fXth[20] = {0};
      //float fXtc[20] = {0};

      unsigned int adcSum = 0.0;
      
      unsigned int maxAdc = 0.0;
      int iphimin = 6666, iphimax = -1;
      int ilaymin = 6666, ilaymax = -1;
      float itmin = 66666666.6, itmax = -6666666666.6;
      
      unsigned int nHits = las_clus->getNhits();
      std::cout << "nhits " << nHits << std::endl;
      if(adcSum>1)
	puts("strange");
      
      for(unsigned int i  = 0; i< nHits;i++){
	float coords[3] = {las_clus->getHitX(i), las_clus->getHitY(i), las_clus->getHitZ(i)};
	int iphi = las_clus->getHitIPhi(i);
	std::cout << " iphi " << iphi << std::endl;

	int it = las_clus->getHitIT(i);
	int ilay = las_clus->getHitLayer(i);
	unsigned int adc = las_clus->getHitAdc(i);
	if(iphi<iphimin){iphimin = iphi;}
	if(iphi>iphimax){iphimax = iphi;}
	if(ilay<ilaymin){ilaymin = ilay;}
	if(ilay>ilaymax){ilaymax = ilay;}
	if(it<itmin){itmin = it;}
	if(it>itmax){itmax = it;}
	adcSum += adc;
	if(adc>maxAdc){
	  maxAdc = adc;
	}
	//	printf("test %d %f %f %f %d %d %d %d %f %f %d %f %f\n", _nevent, coords[0],coords[1],coords[2],adc,ilay,iphi,it,interceptxyff,slopexyff,trk_nclus,zfirst,zlast);

	m_nevent = _nevent;
	m_hitx = coords[0];
	m_hity = coords[1];
	m_hitz = coords[2];
	m_hitadc = adc;
	m_hitlayer = ilay;
	m_hitpad = iphi;
	m_hittbin = it;
	m_interxy = interceptxyff;
	m_slopexy = slopexyff;
	m_interxz = interceptxzff;
	m_slopexz = slopexzff;
	m_interyz = interceptyzff;
	m_slopeyz = slopeyzff;
	m_nclus = trk_nclus;
	m_zfirst = zfirst;
	m_zlast = zlast;
	
	m_hittree->Fill();

	/*
	//"laser hitinfo";"ev:x:y:z:adc:lay:phibin:zbin:intfit:slfit:nclus:zfirst:zlast");
	fXth[nth++] = _nevent;
	
	fXth[nth++] = las_clus->getHitX(i);
	fXth[nth++] = las_clus->getHitY(i);
	fXth[nth++] = las_clus->getHitZ(i);
	
	fXth[nth++] = adc;
	fXth[nth++] = ilay;
	fXth[nth++] = iphi;
	fXth[nth++] = it;
	fXth[nth++] = interceptff;
	fXth[nth++] = slopeff;
	
	fXth[nth++] = trk_nclus;
	fXth[nth++] = zfirst;
	fXth[nth++] = zlast;
	
	_ntp_trk_hit->Fill(_nevent, las_clus->getHitX(i),las_clus->getHitX(i),las_clus->getHitX(i),adc,ilay,iphi,it,interceptff,slopeff,trk_nclus,zfirst,zlast);
		if(_nevent>1000000){
	_ntp_trk_hit->Fill(fXth);
	*/
      }
	
    
    
      int phisize =  iphimax - iphimin + 1;
      int lsize =  ilaymax - ilaymin + 1;
      int tsize = itmax - itmin +1;
      
      m_clux = tcx;
      m_cluy = tcy;
      m_cluz = tcz;
      m_cluadc = adcSum;
      m_clumaxadc = maxAdc;
      m_size = nHits;
      m_sizel = lsize;
      m_sizephi =  phisize;
      m_sizet = tsize;

      m_clustree->Fill();
      
      
      /*
       */
    }
  }
  cout << "number of seeds is " << numberofseeds << endl;
  if(_write_ntp){
    _ntp_max->Fill(_nevent,intmin,intmax,slmin,slmax);
  }

  
  cout << "return " << numberofseeds << endl;
      /*
    if(!keep_event){
    cout << " ABORT !keep " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
    }
  */
  return Fun4AllReturnCodes::EVENT_OK;
}


int PHLineLaserReco::Setup(PHCompositeNode* topNode)
{
  cout << "Called Setup" << endl;
  cout << "topNode:" << topNode << endl;
  // PHTrackSeeding::Setup(topNode);
  //  int ret = GetNodes(topNode);
  //return ret;
  GetNodes(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHLineLaserReco::End(PHCompositeNode* topNode)
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
    _ntp_trk->Write();
    m_hittree->Write();
    m_clustree->Write();
    _tfile->Close();
    cout << "Called End " << endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}


