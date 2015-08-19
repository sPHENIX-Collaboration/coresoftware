#include "PHG4SvtxTrackProjection.h"
#include "SvtxTrackMap.h"
#include "SvtxTrack.h"
#include "PHG4HoughTransform.h"

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>

// PHENIX Geant4 includes
#include <g4cemc/RawTowerGeom.h>
#include <g4cemc/RawTowerContainer.h>
#include <g4cemc/RawTower.h>
#include <g4cemc/RawClusterContainer.h>
#include <g4cemc/RawCluster.h>

// standard includes
#include <iostream>
#include <vector>

using namespace std;

PHG4SvtxTrackProjection::PHG4SvtxTrackProjection(const string &name) :
  SubsysReco(name),
  _num_cal_layers(4),
  _mag_extent(156.5) // middle of Babar magent
{
  _cal_radii.assign(_num_cal_layers,NAN);
  _cal_names.push_back("PRES"); // PRES not yet in G4
  _cal_names.push_back("CEMC");
  _cal_names.push_back("HCALIN");
  _cal_names.push_back("HCALOUT");
}

int PHG4SvtxTrackProjection::Init(PHCompositeNode *topNode) 
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxTrackProjection::InitRun(PHCompositeNode *topNode) 
{
  for (int i=0; i<_num_cal_layers; ++i) {
    string nodename = "TOWERGEOM_" + _cal_names[i];
    RawTowerGeom *geo = findNode::getClass<RawTowerGeom>(topNode,nodename.c_str());
    if (geo) _cal_radii[i] = geo->get_radius();
  }
  
  if (verbosity >= 0) {
    cout << "================== PHG4SvtxTrackProjection::InitRun() =====================" << endl;
    cout << " CVS Version: $Id: PHG4SvtxTrackProjection.C,v 1.10 2015/04/21 23:47:10 pinkenbu Exp $" << endl;
    for (int i=0;i<_num_cal_layers;++i) {
      if (!isnan(_cal_radii[i])) {
	cout << " " << _cal_names[i] << " projection radius: " << _cal_radii[i] << " cm" << endl;
      }
    }
    cout << " projections still curl after the mag field" << endl;
    cout << " projections start from the vertex momentum vector (M.S. effects could be large)" << endl;
    cout << " projections don't correct for the slat HCAL geometry" << endl;
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxTrackProjection::process_event(PHCompositeNode *topNode)
{
  if(verbosity > 1) cout << "PHG4SvtxTrackProjection::process_event -- entered" << endl;

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  // Pull the reconstructed track information off the node tree...
  SvtxTrackMap* _g4tracks = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_g4tracks) {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  for (int i=0;i<_num_cal_layers;++i) {

    if (isnan(_cal_radii[i])) continue;

    if (verbosity > 1) cout << "Projecting tracks into: " << _cal_names[i] << endl;

    // pull the tower geometry
    string towergeonodename = "TOWERGEOM_" + _cal_names[i];
    RawTowerGeom *towergeo = findNode::getClass<RawTowerGeom>(topNode,towergeonodename.c_str());
    if (!towergeo) {
      cerr << PHWHERE << " ERROR: Can't find node " << towergeonodename << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    // pull the towers
    string towernodename = "TOWER_" + _cal_names[i];
    RawTowerContainer *towerList = findNode::getClass<RawTowerContainer>(topNode,towernodename.c_str());
    if (!towerList) {
      cerr << PHWHERE << " ERROR: Can't find node " << towernodename << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    // pull the clusters
    string clusternodename = "CLUSTER_" + _cal_names[i];
    RawClusterContainer *clusterList = findNode::getClass<RawClusterContainer>(topNode,clusternodename.c_str());
    if (!clusterList) {
      cerr << PHWHERE << " ERROR: Can't find node " << clusternodename << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }    
    
    // loop over all tracks
    for (SvtxTrackMap::Iter iter = _g4tracks->begin();
	 iter != _g4tracks->end();
	 ++iter) {
      SvtxTrack *track = &iter->second;

      if (verbosity > 1) cout << "projecting track id " << track->getTrackID() << endl;

      if (verbosity > 1) {
	cout << " track pt = " << sqrt(pow(track->get3Momentum(0),2) +
				       pow(track->get3Momentum(1),2)) << endl;
      }

      // curved tracks inside mag field
      // straight projections thereafter
      std::vector<double> point;
      point.assign(3,-9999.);
      //if (_cal_radii[i] < _mag_extent) {
      // curved projections inside field
      _hough.projectToRadius(*track,_cal_radii[i],point);

      if (isnan(point[0])) continue;
      if (isnan(point[1])) continue;
      if (isnan(point[2])) continue;
      // } else {
      // 	// straight line projections after mag field exit
      // 	_hough.projectToRadius(track,_mag_extent-0.05,point);
      // 	if (isnan(point[0])) continue;
      // 	if (isnan(point[1])) continue;
      // 	if (isnan(point[2])) continue;

      // 	std::vector<double> point2;
      // 	point2.assign(3,-9999.);
      // 	_hough.projectToRadius(track,_mag_extent+0.05,point2);
      // 	if (isnan(point2[0])) continue;
      // 	if (isnan(point2[1])) continue;
      // 	if (isnan(point2[2])) continue;

      // 	// find intersection of r and z


      // find x,y of intersection
      //}
      double x = point[0];
      double y = point[1];
      double z = point[2];

      double phi = atan2(y,x);
      double eta = asinh(z/sqrt(x*x+y*y));

      if (verbosity > 1) {
	cout << " initial track phi = " << atan2(track->get3Momentum(1),track->get3Momentum(0));
	cout << ", eta = " << asinh(track->get3Momentum(2)/sqrt(pow(track->get3Momentum(0),2)+pow(track->get3Momentum(1),2))) << endl;
	cout << " calorimeter phi = " << phi << ", eta = " << eta << endl;
      }

      // projection is outside the detector extent
      // \todo towergeo doesn't make this easy to extract, but this should be
      // fetched from the node tree instead of hardcoded
      if (fabs(eta) >= 1.1) continue;

      // calculate 3x3 tower energy
      int binphi = towergeo->get_phibin(phi);
      int bineta = towergeo->get_etabin(eta);

      double energy_3x3 = 0.0;
      for (int iphi = binphi-1; iphi < binphi+2; ++iphi) { 
      	for (int ieta = bineta-1; ieta < bineta+2; ++ieta) { 

	  // wrap around
	  int wrapphi = iphi;
	  if (wrapphi < 0) {
	    wrapphi = towergeo->get_phibins() + wrapphi;
	  }
	  if (wrapphi >= towergeo->get_phibins()) {
	    wrapphi = wrapphi - towergeo->get_phibins();
	  }

	  // edges
	  if (ieta < 0) continue;
	  if (ieta >= towergeo->get_etabins()) continue;

	  RawTower* tower = towerList->getTower(ieta,wrapphi);
	  if (tower) {
	    energy_3x3 += tower->get_energy();

	    if (verbosity > 1) cout << " tower " << ieta << " " << wrapphi << " energy = " << tower->get_energy() << endl;
	  }
      	}
      }

      track->set_cal_energy_3x3(i,energy_3x3);

      // loop over all clusters and find nearest
      double min_r = DBL_MAX;
      double min_index = -9999;
      double min_dphi = NAN;
      double min_deta = NAN;
      double min_e = NAN;
      for (unsigned int k = 0; k < clusterList->size(); ++k) {

	RawCluster *cluster = clusterList->getCluster(k);

	double dphi = atan2(sin(phi-cluster->get_phi()),cos(phi-cluster->get_phi()));
	double deta = eta-cluster->get_eta();
	double r = sqrt(pow(dphi,2)+pow(deta,2));

	if (r < min_r) {
	  min_index = k;
	  min_r = r;
	  min_dphi = dphi;
	  min_deta = deta;
	  min_e = cluster->get_energy();
	}
      }

      if (min_index != -9999) {
	track->set_cal_dphi(i,min_dphi);
	track->set_cal_deta(i,min_deta);
	track->set_cal_cluster_id(i,min_index);
	track->set_cal_cluster_e(i,min_e);

	if (verbosity > 1) {
	  cout << " nearest cluster dphi = " << min_dphi << " deta = " << min_deta << " e = " << min_e << endl;
	}
      }

    } // end track loop
  } // end calorimeter layer loop

 
  if(verbosity > 1) cout << "PHG4SvtxTrackProjection::process_event -- exited" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxTrackProjection::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
