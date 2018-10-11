/*!
 *  \file		PHG4GenFitTrackProjection.C
 *  \brief		Projects into calorimeters and fills track cal fields using GenFit
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include "PHG4GenFitTrackProjection.h"
#include "SvtxTrackMap.h"
#include "SvtxTrack.h"
#include "PHG4HoughTransform.h"

#include <phgenfit/Fitter.h>
#include <phgenfit/PlanarMeasurement.h>
#include <phgenfit/Track.h>
#include <phgenfit/SpacepointMeasurement.h>

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phgeom/PHGeomUtility.h>
#include <phfield/PHFieldUtility.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <GenFit/RKTrackRep.h>
#include <GenFit/FieldManager.h>

// PHENIX Geant4 includes
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterUtility.h>

//ROOT
#include <TGeoManager.h>

// standard includes
#include <iostream>
#include <vector>
#include <memory>

//#define DEBUG

#define LogDebug(exp)		std::cout<<"DEBUG: "	<<__FILE__<<": "<<__LINE__<<": "<< exp <<std::endl
#define LogError(exp)		std::cout<<"ERROR: "	<<__FILE__<<": "<<__LINE__<<": "<< exp <<std::endl
#define LogWarning(exp)	std::cout<<"WARNING: "	<<__FILE__<<": "<<__LINE__<<": "<< exp <<std::endl

using namespace std;

PHG4GenFitTrackProjection::PHG4GenFitTrackProjection(const string &name, const int pid_guess) :
		SubsysReco(name),

		_fitter(nullptr),

		_pid_guess(pid_guess),

		_num_cal_layers(4)
{
	_cal_radii.assign(_num_cal_layers, NAN);
	_cal_names.push_back("PRES"); // PRES not yet in G4
	_cal_names.push_back("CEMC");
	_cal_names.push_back("HCALIN");
	_cal_names.push_back("HCALOUT");
	_cal_types.push_back(SvtxTrack::PRES); // PRES not yet in G4
	_cal_types.push_back(SvtxTrack::CEMC);
	_cal_types.push_back(SvtxTrack::HCALIN);
	_cal_types.push_back(SvtxTrack::HCALOUT);
}

int PHG4GenFitTrackProjection::Init(PHCompositeNode *topNode) {
	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4GenFitTrackProjection::InitRun(PHCompositeNode *topNode) {
	for (int i = 0; i < _num_cal_layers; ++i) {
		string nodename = "TOWERGEOM_" + _cal_names[i];
		RawTowerGeomContainer *geo = findNode::getClass<RawTowerGeomContainer>(
				topNode, nodename.c_str());
		if (geo)
			_cal_radii[i] = geo->get_radius();//+0.5*geo->get_thickness();
	}

	TGeoManager* tgeo_manager = PHGeomUtility::GetTGeoManager(topNode);

#ifdef DEBUG
	tgeo_manager->Export("Geo_extract.root");
#endif

  PHField * field = PHFieldUtility::GetFieldMapNode(nullptr, topNode);
	_fitter = PHGenFit::Fitter::getInstance(tgeo_manager,field,
			"DafRef",
			"RKTrackRep", false);

	_fitter->set_verbosity(Verbosity());

	if (!_fitter) {
		cerr << PHWHERE << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}

	if (Verbosity() > 0) {
		cout
				<< "================== PHG4GenFitTrackProjection::InitRun() ====================="
				<< endl;
		for (int i = 0; i < _num_cal_layers; ++i) {
			if (!std::isnan(_cal_radii[i])) {
				cout << " " << _cal_names[i] << " projection radius: "
						<< _cal_radii[i] << " cm" << endl;
			}
		}
		cout << " projections still curl after the mag field" << endl;
		cout
				<< " projections start from the vertex momentum vector (M.S. effects could be large)"
				<< endl;
		cout << " projections don't correct for the slat HCAL geometry" << endl;
		cout
				<< "==========================================================================="
				<< endl;
	}

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4GenFitTrackProjection::process_event(PHCompositeNode *topNode) {
	if (Verbosity() > 1)
		cout << "PHG4GenFitTrackProjection::process_event -- entered" << endl;

	//---------------------------------
	// Get Objects off of the Node Tree
	//---------------------------------

	// Pull the reconstructed track information off the node tree...
	SvtxTrackMap* _g4tracks = findNode::getClass<SvtxTrackMap>(topNode,
			"SvtxTrackMap");
	if (!_g4tracks) {
		cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap." << endl;
		return Fun4AllReturnCodes::ABORTRUN;
	}

	for (int i = 0; i < _num_cal_layers; ++i) {

		if (std::isnan(_cal_radii[i]))
			continue;

		if (Verbosity() > 1)
			cout << "Projecting tracks into: " << _cal_names[i] << endl;

		// pull the tower geometry
		string towergeonodename = "TOWERGEOM_" + _cal_names[i];
		RawTowerGeomContainer *towergeo = findNode::getClass<
				RawTowerGeomContainer>(topNode, towergeonodename.c_str());
		if (!towergeo) {
			cerr << PHWHERE << " ERROR: Can't find node " << towergeonodename
					<< endl;
			return Fun4AllReturnCodes::ABORTRUN;
		}

		// pull the towers
		string towernodename = "TOWER_CALIB_" + _cal_names[i];
		RawTowerContainer *towerList = findNode::getClass<RawTowerContainer>(
				topNode, towernodename.c_str());
		if (!towerList) {
			cerr << PHWHERE << " ERROR: Can't find node " << towernodename
					<< endl;
			return Fun4AllReturnCodes::ABORTRUN;
		}

		// pull the clusters
		string clusternodename = "CLUSTER_" + _cal_names[i];
		RawClusterContainer *clusterList = findNode::getClass<
				RawClusterContainer>(topNode, clusternodename.c_str());
		if (!clusterList) {
			cerr << PHWHERE << " ERROR: Can't find node " << clusternodename
					<< endl;
			return Fun4AllReturnCodes::ABORTRUN;
		}

		// loop over all tracks
		for (SvtxTrackMap::Iter iter = _g4tracks->begin();
				iter != _g4tracks->end(); ++iter) {
			SvtxTrack *track = iter->second;
#ifdef DEBUG
			cout
			<<__LINE__
			<<": track->get_charge(): "<<track->get_charge()
			<<endl;
#endif
			if(!track) {
				if(Verbosity() >= 2) LogWarning("!track");
				continue;
			}

			if (Verbosity() > 1)
				cout << "projecting track id " << track->get_id() << endl;

			if (Verbosity() > 1) {
				cout << " track pt = " << track->get_pt() << endl;
			}

			std::vector<double> point;
			point.assign(3, -9999.);

			auto last_state_iter = --track->end_states();

			SvtxTrackState * trackstate = last_state_iter->second;

			if(!trackstate) {
				if(Verbosity() >= 2) LogWarning("!trackstate");
				continue;
			}

			auto pdg = unique_ptr<TDatabasePDG> (TDatabasePDG::Instance());
			int reco_charge = track->get_charge();
			int gues_charge = pdg->GetParticle(_pid_guess)->Charge();
			if(reco_charge*gues_charge<0) _pid_guess *= -1;
#ifdef DEBUG
			cout
			<<__LINE__
			<<": guess charge: " << gues_charge
			<<": reco charge: " << reco_charge
			<<": pid: " << _pid_guess
			<<": pT: " << sqrt(trackstate->get_px()*trackstate->get_px() + trackstate->get_py()*trackstate->get_py())
			<<endl;
#endif

			auto rep = unique_ptr<genfit::AbsTrackRep> (new genfit::RKTrackRep(_pid_guess));

			unique_ptr<genfit::MeasuredStateOnPlane> msop80 = nullptr;

			{
				TVector3 pos(trackstate->get_x(), trackstate->get_y(), trackstate->get_z());
				//pos.SetXYZ(0.01,0,0);

				TVector3 mom(trackstate->get_px(), trackstate->get_py(), trackstate->get_pz());
				//mom.SetXYZ(1,0,0);

				TMatrixDSym cov(6);
				for (int i = 0; i < 6; ++i) {
					for (int j = 0; j < 6; ++j) {
						cov[i][j] = trackstate->get_error(i, j);
					}
				}

				msop80 = unique_ptr<genfit::MeasuredStateOnPlane> (new genfit::MeasuredStateOnPlane(rep.get()));

				msop80->setPosMomCov(pos, mom, cov);
			}

#ifdef DEBUG
			{
				double x = msop80->getPos().X();
				double y = msop80->getPos().Y();
				double z = msop80->getPos().Z();
//				double px = msop80->getMom().X();
//				double py = msop80->getMom().Y();
				double pz = msop80->getMom().Z();
				genfit::FieldManager *field_mgr = genfit::FieldManager::getInstance();
				double Bx=0, By=0, Bz=0;
				field_mgr->getFieldVal(x,y,z,Bx,By,Bz);
				cout
				<< __LINE__
				<< ": { " << msop80->getPos().Perp() << ", " << msop80->getPos().Phi() << ", " << msop80->getPos().Eta() << "} @ "
				//<< "{ " << Bx << ", " << By << ", " << Bz << "}"
				<< "{ " << msop80->getMom().Perp() << ", " << msop80->getMom().Phi() << ", " << pz << "} "
				<<endl;
				//msop80->Print();
			}
#endif
			try {
				rep->extrapolateToCylinder(*msop80, _cal_radii[i], TVector3(0,0,0),  TVector3(0,0,1));
				//rep->extrapolateToCylinder(*msop80, 5., TVector3(0,0,0),  TVector3(0,0,1));
			} catch (...) {
				if(Verbosity() >= 2) LogWarning("extrapolateToCylinder failed");
				continue;
			}

#ifdef DEBUG
			{
				cout<<__LINE__<<endl;
				//msop80->Print();
				double x = msop80->getPos().X();
				double y = msop80->getPos().Y();
				double z = msop80->getPos().Z();
//				double px = msop80->getMom().X();
//				double py = msop80->getMom().Y();
				double pz = msop80->getMom().Z();
				genfit::FieldManager *field_mgr = genfit::FieldManager::getInstance();
				double Bx=0, By=0, Bz=0;
				field_mgr->getFieldVal(x,y,z,Bx,By,Bz);
				cout
				<< __LINE__
				<< ": { " << msop80->getPos().Perp() << ", " << msop80->getPos().Phi() << ", " << msop80->getPos().Eta() << "} @ "
				//<< "{ " << Bx << ", " << By << ", " << Bz << "}"
				<< "{ " << msop80->getMom().Perp() << ", " << msop80->getMom().Phi() << ", " << pz << "} "
				<<endl;
			}
#endif

			point[0] = msop80->getPos().X();
			point[1] = msop80->getPos().Y();
			point[2] = msop80->getPos().Z();

#ifdef DEBUG
			cout
			<<__LINE__
			<<": GenFit: {"
			<< point[0] <<", "
			<< point[1] <<", "
			<< point[2] <<" }"
			<<endl;
#endif

			if (std::isnan(point[0]))
				continue;
			if (std::isnan(point[1]))
				continue;
			if (std::isnan(point[2]))
				continue;

			double x = point[0];
			double y = point[1];
			double z = point[2];

			double phi = atan2(y, x);
			double eta = asinh(z / sqrt(x * x + y * y));

			if (Verbosity() > 1) {
				cout << " initial track phi = " << track->get_phi();
				cout << ", eta = " << track->get_eta() << endl;
				cout << " calorimeter phi = " << phi << ", eta = " << eta
						<< endl;
			}

			// projection is outside the detector extent
			// TODO towergeo doesn't make this easy to extract, but this should be
			// fetched from the node tree instead of hardcoded
			if (fabs(eta) >= 1.0)
				continue;

			// calculate 3x3 tower energy
			int binphi = towergeo->get_phibin(phi);
			int bineta = towergeo->get_etabin(eta);

			double energy_3x3 = 0.0;
			double energy_5x5 = 0.0;
			for (int iphi = binphi - 2; iphi <= binphi + 2; ++iphi) {
				for (int ieta = bineta - 2; ieta <= bineta + 2; ++ieta) {

					// wrap around
					int wrapphi = iphi;
					if (wrapphi < 0) {
						wrapphi = towergeo->get_phibins() + wrapphi;
					}
					if (wrapphi >= towergeo->get_phibins()) {
						wrapphi = wrapphi - towergeo->get_phibins();
					}

					// edges
					if (ieta < 0)
						continue;
					if (ieta >= towergeo->get_etabins())
						continue;

					RawTower* tower = towerList->getTower(ieta, wrapphi);
					if (tower) {

						energy_5x5 += tower->get_energy();
						if (abs(iphi - binphi) <= 1 and abs(ieta - bineta) <= 1)
							energy_3x3 += tower->get_energy();

						if (Verbosity() > 1)
							cout << " tower " << ieta << " " << wrapphi
									<< " energy = " << tower->get_energy()
									<< endl;
					}
				}
			}

			track->set_cal_energy_3x3(_cal_types[i], energy_3x3);
			track->set_cal_energy_5x5(_cal_types[i], energy_5x5);

			// loop over all clusters and find nearest
			double min_r = DBL_MAX;
			double min_index = -9999;
			double min_dphi = NAN;
			double min_deta = NAN;
			double min_e = NAN;
#ifdef DEBUG
			double min_cluster_phi = NAN;
#endif
			for (const auto & iterator : clusterList->getClustersMap()) {

				const RawCluster *cluster = iterator.second;

				//! eta as location mark of cluster relative to (0,0,0)
				const float cluster_eta = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0,0,0));

				double dphi = atan2(sin(phi - cluster->get_phi()),
						cos(phi - cluster->get_phi()));
				double deta = eta - cluster_eta;
				double r = sqrt(pow(dphi, 2) + pow(deta, 2));

				if (r < min_r) {
					min_index = iterator.first;
					min_r = r;
					min_dphi = dphi;
					min_deta = deta;
					min_e = cluster->get_energy();
#ifdef DEBUG
					min_cluster_phi = cluster->get_phi();
#endif
				}
			}

			if (min_index != -9999) {
				track->set_cal_dphi(_cal_types[i], min_dphi);
				track->set_cal_deta(_cal_types[i], min_deta);
				track->set_cal_cluster_id(_cal_types[i], min_index);
				track->set_cal_cluster_e(_cal_types[i], min_e);

#ifdef DEBUG
			cout
			<<__LINE__
			<<": min_cluster_phi: "<<min_cluster_phi
			<<endl;
#endif

				if (Verbosity() > 1) {
					cout << " nearest cluster dphi = " << min_dphi << " deta = "
							<< min_deta << " e = " << min_e << endl;
				}
			}

		} // end track loop
	} // end calorimeter layer loop

	if (Verbosity() > 1)
		cout << "PHG4GenFitTrackProjection::process_event -- exited" << endl;

	return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4GenFitTrackProjection::End(PHCompositeNode *topNode) {
	return Fun4AllReturnCodes::EVENT_OK;
}
