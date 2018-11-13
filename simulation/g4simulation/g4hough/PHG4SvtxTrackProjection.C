#include "PHG4SvtxTrackProjection.h"
#include "SvtxTrackMap.h"
#include "SvtxTrack.h"
#include "PHG4HoughTransform.h"

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

// PHENIX Geant4 includes
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterUtility.h>

// standard includes
#include <iostream>
#include <vector>

using namespace std;

PHG4SvtxTrackProjection::PHG4SvtxTrackProjection(const string &name) :
  SubsysReco(name),
  _num_cal_layers(4),
  _magfield(1.5),
  _mag_extent(156.5) // middle of Babar magent
{
  _cal_radii.assign(_num_cal_layers,NAN);
  _cal_names.push_back("PRES"); // PRES not yet in G4
  _cal_names.push_back("CEMC");
  _cal_names.push_back("HCALIN");
  _cal_names.push_back("HCALOUT");
  _cal_types.push_back(SvtxTrack::PRES); // PRES not yet in G4
  _cal_types.push_back(SvtxTrack::CEMC);
  _cal_types.push_back(SvtxTrack::HCALIN);
  _cal_types.push_back(SvtxTrack::HCALOUT);
}

int PHG4SvtxTrackProjection::Init(PHCompositeNode *topNode) 
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxTrackProjection::InitRun(PHCompositeNode *topNode) 
{
  for (int i=0; i<_num_cal_layers; ++i) {
    string nodename = "TOWERGEOM_" + _cal_names[i];
    RawTowerGeomContainer *geo = findNode::getClass<RawTowerGeomContainer>(topNode,nodename.c_str());
    if (geo) _cal_radii[i] = geo->get_radius();
  }
  
  if (Verbosity() > 0) {
    cout << "================== PHG4SvtxTrackProjection::InitRun() =====================" << endl;
    for (int i=0;i<_num_cal_layers;++i) {
      if (!std::isnan(_cal_radii[i])) {
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
  if(Verbosity() > 1) cout << "PHG4SvtxTrackProjection::process_event -- entered" << endl;

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

    if (std::isnan(_cal_radii[i])) continue;

    if (Verbosity() > 1) cout << "Projecting tracks into: " << _cal_names[i] << endl;

    // pull the tower geometry
    string towergeonodename = "TOWERGEOM_" + _cal_names[i];
    RawTowerGeomContainer *towergeo = findNode::getClass<RawTowerGeomContainer>(topNode,towergeonodename.c_str());
    if (!towergeo) {
      cerr << PHWHERE << " ERROR: Can't find node " << towergeonodename << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    // pull the towers
    string towernodename = "TOWER_CALIB_" + _cal_names[i];
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
      SvtxTrack *track = iter->second;

      if (Verbosity() > 1) cout << "projecting track id " << track->get_id() << endl;

      if (Verbosity() > 1) {
	cout << " track pt = " << track->get_pt() << endl;
      }

      // curved tracks inside mag field
      // straight projections thereafter
      std::vector<double> point;
      point.assign(3,-9999.);
      //if (_cal_radii[i] < _mag_extent) {
      // curved projections inside field

      _hough.projectToRadius(track,_magfield,_cal_radii[i],point);

      if (std::isnan(point[0])) continue;
      if (std::isnan(point[1])) continue;
      if (std::isnan(point[2])) continue;
      // } else {
      // 	// straight line projections after mag field exit
      // 	_hough.projectToRadius(track,_mag_extent-0.05,point);
      // 	if (std::isnan(point[0])) continue;
      // 	if (std::isnan(point[1])) continue;
      // 	if (std::isnan(point[2])) continue;

      // 	std::vector<double> point2;
      // 	point2.assign(3,-9999.);
      // 	_hough.projectToRadius(track,_mag_extent+0.05,point2);
      // 	if (std::isnan(point2[0])) continue;
      // 	if (std::isnan(point2[1])) continue;
      // 	if (std::isnan(point2[2])) continue;

      // 	// find intersection of r and z


      // find x,y of intersection
      //}
      double x = point[0];
      double y = point[1];
      double z = point[2];

      double phi = atan2(y,x);
      double eta = asinh(z/sqrt(x*x+y*y));

      if (Verbosity() > 1) {
	cout << " initial track phi = " << track->get_phi();
	cout << ", eta = " << track->get_eta() << endl;
	cout << " calorimeter phi = " << phi << ", eta = " << eta << endl;
      }

      // projection is outside the detector extent
      // \todo towergeo doesn't make this easy to extract, but this should be
      // fetched from the node tree instead of hardcoded
      if (fabs(eta) >= 1.0) continue;

      // calculate 3x3 tower energy
      int binphi = towergeo->get_phibin(phi);
      int bineta = towergeo->get_etabin(eta);

      double energy_3x3 = 0.0;
      double energy_5x5 = 0.0;
      for (int iphi = binphi-2; iphi <= binphi+2; ++iphi) {
      	for (int ieta = bineta-2; ieta <= bineta+2; ++ieta) {

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

        energy_5x5 += tower->get_energy();
	      if (abs(iphi - binphi)<=1 and abs(ieta - bineta)<=1 )
	        energy_3x3 += tower->get_energy();

	    if (Verbosity() > 1) cout << " tower " << ieta << " " << wrapphi << " energy = " << tower->get_energy() << endl;
	  }
      	}
      }

      track->set_cal_energy_3x3(_cal_types[i],energy_3x3);
      track->set_cal_energy_5x5(_cal_types[i],energy_5x5);

      // loop over all clusters and find nearest
      double min_r = DBL_MAX;
      double min_index = -9999;
      double min_dphi = NAN;
      double min_deta = NAN;
      double min_e = NAN;
      for (const auto & iterator : clusterList->getClustersMap()) {

        const RawCluster *cluster = iterator.second;

  //! eta as location mark of cluster relative to (0,0,0)
  const float cluster_eta = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0,0,0));

	double dphi = atan2(sin(phi-cluster->get_phi()),cos(phi-cluster->get_phi()));
	double deta = eta-cluster_eta;
	double r = sqrt(pow(dphi,2)+pow(deta,2));

	if (r < min_r) {
	  min_index = iterator.first;
	  min_r = r;
	  min_dphi = dphi;
	  min_deta = deta;
	  min_e = cluster->get_energy();
	}
      }

      if (min_index != -9999) {
	track->set_cal_dphi(_cal_types[i],min_dphi);
	track->set_cal_deta(_cal_types[i],min_deta);
	track->set_cal_cluster_id(_cal_types[i],min_index);
	track->set_cal_cluster_e(_cal_types[i],min_e);

	if (Verbosity() > 1) {
	  cout << " nearest cluster dphi = " << min_dphi << " deta = " << min_deta << " e = " << min_e << endl;
	}
      }

    } // end track loop
  } // end calorimeter layer loop

 
  if(Verbosity() > 1) cout << "PHG4SvtxTrackProjection::process_event -- exited" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxTrackProjection::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
