/*!
 *  \file		Track.cc
 *  \brief		Data structure and output of the fitting.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

//STL

//BOOST
//#include <boost/foreach.hpp>

//GenFit
#include <GenFit/Track.h>
#include <GenFit/MeasurementOnPlane.h>
#include <GenFit/KalmanFittedStateOnPlane.h>
#include <GenFit/KalmanFitterInfo.h>
#include <GenFit/KalmanFitter.h>
#include <GenFit/Tools.h>
#include <GenFit/AbsHMatrix.h>

//PHGenFit
#include "Track.h"
#include "Measurement.h"

#define LogDebug(exp)		std::cout<<"DEBUG: "	<<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogError(exp)		std::cout<<"ERROR: "	<<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogWarning(exp)	std::cout<<"WARNING: "	<<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl

#define WILD_DOUBLE -999999

#define _DEBUG_
//#define _PRINT_MATRIX_

namespace PHGenFit {

Track::Track(genfit::AbsTrackRep *rep, TVector3 seed_pos, TVector3 seed_mom, TMatrixDSym seed_cov)
{
	//TODO Add input param check

	verbosity = 2;

	genfit::MeasuredStateOnPlane seedMSoP(rep);
	seedMSoP.setPosMomCov(seed_pos, seed_mom, seed_cov);
	const genfit::StateOnPlane seedSoP(seedMSoP);

	TVectorD seedState(6);
	TMatrixDSym seedCov(6);
	seedMSoP.get6DStateCov(seedState, seedCov);


	_track = new genfit::Track(rep, seedState, seedCov);
	//_track = NEW(genfit::Track)(rep, seedState, seedCov);
}

Track::Track(const PHGenFit::Track &t) {
	_track = new genfit::Track(*(t.getGenFitTrack()));
	_clusterIDs = t.get_cluster_IDs();
}

int Track::addMeasurement(PHGenFit::Measurement*measurement) {

	std::vector<genfit::AbsMeasurement*> msmts;
	msmts.push_back(measurement->getMeasurement());
	_track->insertPoint(new genfit::TrackPoint(msmts, _track));

	_clusterIDs.push_back(measurement->get_cluster_ID());

	return 0;
}

int Track::addMeasurements(std::vector<PHGenFit::Measurement*> &measurements)
{
	for(PHGenFit::Measurement* measurement : measurements)
	{
		std::vector<genfit::AbsMeasurement*> msmts;
		msmts.push_back(measurement->getMeasurement());
		_track->insertPoint(
				new genfit::TrackPoint(msmts, _track));

		//_measurements.push_back(measurement);
		_clusterIDs.push_back(measurement->get_cluster_ID());
	}

	//measurements.clear();

	return 0;
}

int Track::deleteLastMeasurement() {

	_track->deletePoint(-1);

	_clusterIDs.pop_back();

	return 0;
}

Track::~Track()
{
	delete _track;

//	for(PHGenFit::Measurement* measurement : _measurements)
//	{
//		delete measurement;
//	}
//	_measurements.clear();

	_clusterIDs.clear();
}

double Track::extrapolateToPlane(genfit::MeasuredStateOnPlane& state, TVector3 O, TVector3 n, const int tr_point_id) const
{
	double pathlenth = WILD_DOUBLE;

	genfit::SharedPlanePtr destPlane(new genfit::DetPlane(O, n));

	genfit::AbsTrackRep* rep = _track->getCardinalRep();
	genfit::TrackPoint* tp = _track->getPointWithMeasurementAndFitterInfo(
			tr_point_id, rep);
	if (tp == NULL) {
		std::cout << "Track has no TrackPoint with fitterInfo! \n";
		return WILD_DOUBLE;
	}
	std::unique_ptr<genfit::KalmanFittedStateOnPlane> kfsop (new genfit::KalmanFittedStateOnPlane(
			*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate())));
	// extrapolate back to reference plane.
	try {
		pathlenth = rep->extrapolateToPlane(*kfsop, destPlane);
	} catch (genfit::Exception& e) {
		std::cerr << "Exception, next track" << std::endl;
		std::cerr << e.what();
		//delete kfsop;
		return WILD_DOUBLE;
	}

	state = *dynamic_cast<genfit::MeasuredStateOnPlane*> (kfsop.get());

	//delete kfsop;

	return pathlenth;
}

genfit::MeasuredStateOnPlane* Track::extrapolateToPlane(TVector3 O, TVector3 n, const int tr_point_id) const
{
	genfit::MeasuredStateOnPlane* state = new genfit::MeasuredStateOnPlane();
	double pathlenth = this->extrapolateToPlane(*state, O, n, tr_point_id);
	if(pathlenth <= WILD_DOUBLE) {
		delete state;
		return NULL;
	}
	else
		return state;
}

double Track::extrapolateToLine(genfit::MeasuredStateOnPlane& state, TVector3 line_point, TVector3 line_direction, const int tr_point_id) const
{
	double pathlenth = WILD_DOUBLE;

	genfit::AbsTrackRep* rep = _track->getCardinalRep();
	genfit::TrackPoint* tp = _track->getPointWithMeasurementAndFitterInfo(
			tr_point_id, rep);
	if (tp == NULL) {
		std::cout << "Track has no TrackPoint with fitterInfo! \n";
		return WILD_DOUBLE;
	}
	std::unique_ptr<genfit::KalmanFittedStateOnPlane> kfsop (new genfit::KalmanFittedStateOnPlane(
			*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate())));
	// extrapolate back to reference plane.
	try {
		pathlenth = rep->extrapolateToLine(*kfsop, line_point, line_direction);
	} catch (genfit::Exception& e) {
		std::cerr << "Exception, next track" << std::endl;
		std::cerr << e.what();
		//delete kfsop;
		return WILD_DOUBLE;
	}

	state = *dynamic_cast<genfit::MeasuredStateOnPlane*> (kfsop.get());

	//delete kfsop;

	return pathlenth;
}

genfit::MeasuredStateOnPlane* Track::extrapolateToLine(TVector3 line_point, TVector3 line_direction, const int tr_point_id) const
{
	genfit::MeasuredStateOnPlane* state = new genfit::MeasuredStateOnPlane();
	double pathlenth = this->extrapolateToLine(*state, line_point, line_direction, tr_point_id);
	if(pathlenth <= WILD_DOUBLE) {
		delete state;
		return NULL;
	}
	else
		return state;
}

double Track::extrapolateToCylinder(genfit::MeasuredStateOnPlane& state, double radius, TVector3 line_point, TVector3 line_direction, const int tr_point_id, const int direction) const
{
	assert(direction == 1 or direction == -1);

	double pathlenth = WILD_DOUBLE;

	genfit::AbsTrackRep* rep = _track->getCardinalRep();

//	genfit::TrackPoint* tp = _track->getPointWithMeasurementAndFitterInfo(
//			tr_point_id, rep);
//	if (tp == NULL) {
//		std::cout << "Track has no TrackPoint with fitterInfo! \n";
//		return WILD_DOUBLE;
//	}

	bool have_tp_with_fit_info = false;
	std::unique_ptr<genfit::MeasuredStateOnPlane> kfsop = NULL;
	if (_track->getNumPointsWithMeasurement() > 0) {
#ifdef _DEBUG_
//		std::cout<<__LINE__ <<std::endl;
#endif
		genfit::TrackPoint* tp = _track->getPointWithMeasurement(tr_point_id);
		if (tp == NULL) {
			LogError("tp == NULL!");
			return WILD_DOUBLE;
		}
		if (dynamic_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))) {
			if (direction == 1) {
				if (static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(
						rep))->getForwardUpdate()) {
					have_tp_with_fit_info = true;
					kfsop =
							std::unique_ptr < genfit::MeasuredStateOnPlane
									> (new genfit::KalmanFittedStateOnPlane(
											*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(
													rep))->getForwardUpdate())));
				}
			} else {
				if (static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(
						rep))->getBackwardUpdate()) {
					have_tp_with_fit_info = true;
					kfsop =
							std::unique_ptr < genfit::MeasuredStateOnPlane
									> (new genfit::KalmanFittedStateOnPlane(
											*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(
													rep))->getBackwardUpdate())));
				}
			}

		}
	}

	if (!have_tp_with_fit_info) {
#ifdef _DEBUG_
		std::cout<<__LINE__<< ": !have_tp_with_fit_info" <<std::endl;
#endif
		kfsop = std::unique_ptr < genfit::MeasuredStateOnPlane
				> (new genfit::MeasuredStateOnPlane(rep));
		rep->setPosMomCov(*kfsop, _track->getStateSeed(), _track->getCovSeed());
	}

	if(!kfsop) return pathlenth;
	// extrapolate back to reference plane.
	try {
		//rep->extrapolateToLine(*kfsop, line_point, line_direction);
		pathlenth = rep->extrapolateToCylinder(*kfsop, radius, line_point, line_direction);
	} catch (genfit::Exception& e) {
		if(verbosity > 1) {
			LogWarning("Can't extrapolate to cylinder!");
			std::cerr << e.what();
		}
		return WILD_DOUBLE;
	}

	state = *dynamic_cast<genfit::MeasuredStateOnPlane*> (kfsop.get());

	//delete kfsop;

	return pathlenth;
}

genfit::MeasuredStateOnPlane*  Track::extrapolateToCylinder(double radius, TVector3 line_point, TVector3 line_direction, const int tr_point_id, const int direction) const
{
	assert(direction == 1 or direction == -1);
	genfit::MeasuredStateOnPlane* state = new genfit::MeasuredStateOnPlane();
	double pathlenth = this->extrapolateToCylinder(*state, radius, line_point, line_direction, tr_point_id, direction);
	if(pathlenth <= WILD_DOUBLE) {
		delete state;
		return NULL;
	}
	else
		return state;
}

int Track::updateOneMeasurementKalman(
		const std::vector<PHGenFit::Measurement*>& measurements,
		std::map<double, PHGenFit::Track*>& incr_chi2s_new_tracks,
		const int base_tp_idx,
		const int direction,
		const float blowup_factor,
		const bool use_fitted_state) const {

#ifdef _DEBUG_
	std::cout
	<< __LINE__
	<<" : base_tp_idx: " << base_tp_idx
	<<" : direction: " << direction
	<<" : blowup_factor: " << blowup_factor
	<<" : use_fitted_state: " << use_fitted_state
	<< std::endl;
#endif

	if(measurements.size()==0) return -1;

	for (PHGenFit::Measurement* measurement : measurements) {

		PHGenFit::Track* new_track = NULL;

		new_track = new PHGenFit::Track(*this);

//		if(incr_chi2s_new_tracks.size() == 0)
//			new_track = const_cast<PHGenFit::Track*>(this);
//		else
//			new_track = new PHGenFit::Track(*this);

		genfit::Track *track = new_track->getGenFitTrack();
		genfit::AbsTrackRep* rep = track->getCardinalRep();

		bool newFi(true);
		genfit::TrackPoint *tp_base = NULL;
		std::unique_ptr<genfit::MeasuredStateOnPlane> currentState = NULL;
		genfit::SharedPlanePtr plane = NULL;

		if(track->getNumPointsWithMeasurement() > 0) {
			tp_base = track->getPointWithMeasurement(base_tp_idx);
			newFi = !(tp_base->hasFitterInfo(rep));
			//tp_base->Print();
		}
#ifdef _DEBUG_
		std::cout << __LINE__ << ": " <<"newFi: "<<newFi<<std::endl;
#endif
		if (newFi) {
			currentState = std::unique_ptr < genfit::MeasuredStateOnPlane
					> (new genfit::MeasuredStateOnPlane(rep));
			rep->setPosMomCov(*currentState, track->getStateSeed(),
					track->getCovSeed());
		} else {
			try {
				genfit::KalmanFitterInfo* kfi = static_cast<genfit::KalmanFitterInfo*>(tp_base->getFitterInfo(rep));
				if(!kfi) {
#ifdef _DEBUG_
					LogDebug("!kfi");
#endif
					continue;
				}
//#ifdef _DEBUG_
//				std::cout << __LINE__ << "\n ###################################################################"<<std::endl;
//				kfi->Print();
//				std::cout << __LINE__ << "\n ###################################################################"<<std::endl;
//#endif
				if(use_fitted_state) {
					const genfit::MeasuredStateOnPlane* tempFS = &(kfi->getFittedState(true));
					if(!tempFS) {
#ifdef _DEBUG_
						LogDebug("!tempFS");
#endif
						continue;
					}
					currentState = std::unique_ptr < genfit::MeasuredStateOnPlane> (new genfit::MeasuredStateOnPlane(*tempFS));
				} else {
					genfit::MeasuredStateOnPlane* tempUpdate =  kfi->getUpdate(direction);
					if(!tempUpdate) {
#ifdef _DEBUG_
						LogDebug("!tempUpdate");
#endif
						continue;
					}
					currentState = std::unique_ptr < genfit::MeasuredStateOnPlane> (new genfit::MeasuredStateOnPlane(*tempUpdate));
				}

#ifdef _DEBUG_
//				std::cout << __LINE__ << "\n ###################################################################"<<std::endl;
//				kfi->Print();
//				std::cout << __LINE__ << "\n ###################################################################"<<std::endl;
//				tempFS->Print();
//				std::cout << __LINE__ << "\n ###################################################################"<<std::endl;
//				tempUpdate->Print();
//				std::cout << __LINE__ << "\n ###################################################################"<<std::endl;
#endif

				if(blowup_factor > 1) {
					currentState->blowUpCov(blowup_factor, true, 1e6);
				}

			} catch (genfit::Exception e) {
#ifdef _DEBUG_
		std::cout
		<< __LINE__
		<< ": Fitted state not found!"
		<< std::endl;
		std::cerr<<e.what()<<std::endl;
#endif
				currentState = std::unique_ptr < genfit::MeasuredStateOnPlane
						> (new genfit::MeasuredStateOnPlane(rep));
				rep->setPosMomCov(*currentState, track->getStateSeed(),
						track->getCovSeed());
			}
		}
#ifdef _DEBUG_
		std::cout<< __LINE__ <<std::endl;
#endif
		//std::vector<genfit::AbsMeasurement*> msmts;
		//msmts.push_back(measurement->getMeasurement());

		//genfit::TrackPoint *tp = new genfit::TrackPoint(msmts, track);
		//track->insertPoint(tp); // genfit

		/*!
		 * A new TrackPoint created in addMeasurement
		 * PHGenFit: clusterID also registerd
		 */
		new_track->addMeasurement(measurement);

#ifdef _DEBUG_
		std::cout<< __LINE__ <<": clusterIDs size: "<< new_track->get_cluster_IDs().size() <<std::endl;
#endif

		//! Get the pointer of the TrackPoint just created
		genfit::TrackPoint *tp = new_track->getGenFitTrack()->getPoint(-1);
#ifdef _DEBUG_
		std::cout<< __LINE__ <<std::endl;
#endif
		genfit::KalmanFitterInfo* fi = new genfit::KalmanFitterInfo(tp, rep);
		tp->setFitterInfo(fi);
#ifdef _DEBUG_
		std::cout
		<< __LINE__
		<<": track->getPointWithMeasurement(): " << track->getPointWithMeasurement(-1)
		<< std::endl;
#endif
//		if (track->getNumPointsWithMeasurement() > 0) {
//			tp_base = track->getPointWithMeasurement(-1);
//			if (tp_base->hasFitterInfo(rep)) {
//				std::cout << "TP has FI!" << std::endl;
//			}
//		}
#ifdef _DEBUG_
		std::cout<< __LINE__ <<std::endl;
#endif
		const std::vector<genfit::AbsMeasurement*>& rawMeasurements =
				tp->getRawMeasurements();
		// construct plane with first measurement
		plane = rawMeasurements[0]->constructPlane(*currentState);

		//double extLen = rep->extrapolateToPlane(*state, plane);

		try {
			rep->extrapolateToPlane(*currentState, plane);
		} catch (...) {
			if(verbosity > 1) {
				LogWarning("Can not extrapolate to measuremnt: ") << measurement->get_cluster_ID() <<std::endl;
			}
			continue;
		}
#ifdef _DEBUG_
		std::cout<< __LINE__ <<std::endl;
#endif
		fi->setPrediction(currentState->clone(), direction);
		genfit::MeasuredStateOnPlane *state = fi->getPrediction(direction);
#ifdef _DEBUG_
		std::cout<< __LINE__ <<std::endl;
#endif
		TVectorD stateVector(state->getState());
		TMatrixDSym cov(state->getCov());
#ifdef _DEBUG_
		std::cout<< __LINE__ <<std::endl;
#endif
		for (std::vector<genfit::AbsMeasurement*>::const_iterator it =
				rawMeasurements.begin(); it != rawMeasurements.end(); ++it) {
			fi->addMeasurementsOnPlane(
					(*it)->constructMeasurementsOnPlane(*state));
		}

		double chi2inc = 0;
		double ndfInc = 0;
#ifdef _DEBUG_
		std::cout<< __LINE__ <<std::endl;
#endif
		// update(s)
		const std::vector<genfit::MeasurementOnPlane *>& measurements_on_plane = fi->getMeasurementsOnPlane();
#ifdef _DEBUG_
		std::cout
		<< __LINE__
		<< ": size of fi's MeasurementsOnPlane: " << measurements_on_plane.size()
		<<std::endl;
#endif
		for (std::vector<genfit::MeasurementOnPlane *>::const_iterator it =
				measurements_on_plane.begin();
				it != measurements_on_plane.end(); ++it) {
			const genfit::MeasurementOnPlane& mOnPlane = **it;
			//const double weight = mOnPlane.getWeight();

			const TVectorD& measurement(mOnPlane.getState());
			const genfit::AbsHMatrix* H(mOnPlane.getHMatrix());
			// (weighted) cov
			const TMatrixDSym& V(mOnPlane.getCov()); //Covariance of measurement noise v_{k}

			TVectorD res(measurement - H->Hv(stateVector));

			// If hit, do Kalman algebra.
			{
				// calculate kalman gain ------------------------------
				// calculate covsum (V + HCH^T) and invert
				TMatrixDSym covSumInv(cov);
				H->HMHt(covSumInv);
				covSumInv += V;
				try{
					genfit::tools::invertMatrix(covSumInv);
				} catch (genfit::Exception e) {
#ifdef _DEBUG_
					LogDebug("cannot invert matrix.");
#endif
					continue;
				}

				TMatrixD CHt(H->MHt(cov));
#ifdef _PRINT_MATRIX_
				std::cout <<__LINE__ <<": V_{k}:" << std::endl;
				V.Print();
				std::cout <<__LINE__ <<": R_{k}^{-1}:" << std::endl;
				covSumInv.Print();
				std::cout <<__LINE__ <<": C_{k|k-1}:" << std::endl;
				cov.Print();
				std::cout <<__LINE__ <<": C_{k|k-1} H_{k}^{T} :" << std::endl;
				CHt.Print();
				std::cout <<__LINE__ <<": K_{k} :" << std::endl;
				TMatrixD Kk(CHt, TMatrixD::kMult, covSumInv);
				Kk.Print();
				std::cout <<__LINE__ <<": res:" << std::endl;
				res.Print();
#endif
				TVectorD update(
						TMatrixD(CHt, TMatrixD::kMult, covSumInv) * res);
				//TMatrixD(CHt, TMatrixD::kMult, covSumInv).Print();

				stateVector += update; // x_{k|k} = x_{k|k-1} + K_{k} r_{k|k-1}
				covSumInv.Similarity(CHt); // with (C H^T)^T = H C^T = H C  (C is symmetric)
				cov -= covSumInv; //C_{k|k}
			}

			TVectorD resNew(measurement - H->Hv(stateVector));

			// Calculate chi2
			TMatrixDSym HCHt(cov); //C_{k|k}
			H->HMHt(HCHt);
			HCHt -= V;
			HCHt *= -1;

			try{
				genfit::tools::invertMatrix(HCHt);
			} catch (genfit::Exception e) {
#ifdef _DEBUG_
				LogDebug("cannot invert matrix.");
#endif
				continue;
			}
			chi2inc += HCHt.Similarity(resNew);

			ndfInc += measurement.GetNrows();
#ifdef _PRINT_MATRIX_
				std::cout <<__LINE__ <<": V - HCHt:" << std::endl;
				HCHt.Print();
				std::cout <<__LINE__ <<": resNew:" << std::endl;
				resNew.Print();
#endif

#ifdef _DEBUG_
				std::cout <<__LINE__ <<": ndfInc:  " << ndfInc << std::endl;
				std::cout <<__LINE__ <<": chi2inc: " << chi2inc << std::endl;
#endif
			genfit::KalmanFittedStateOnPlane* updatedSOP = new genfit::KalmanFittedStateOnPlane(
					*state, chi2inc, ndfInc);
			fi->setUpdate(updatedSOP, direction);
		} //loop measurements_on_plane

		//FIXME why chi2 could be smaller than 0?
		if (chi2inc > 0)
			incr_chi2s_new_tracks.insert(std::make_pair(chi2inc,new_track));

	}//loop measurments

	return 0;
}

double Track::extrapolateToPoint(genfit::MeasuredStateOnPlane& state, TVector3 P, const int tr_point_id) const
{
	double pathlenth = WILD_DOUBLE;
	genfit::AbsTrackRep* rep = _track->getCardinalRep();
	genfit::TrackPoint* tp = _track->getPointWithMeasurementAndFitterInfo(
			tr_point_id, rep);
	if (tp == NULL) {
		std::cout << "Track has no TrackPoint with fitterInfo! \n";
		return WILD_DOUBLE;
	}
	std::unique_ptr<genfit::KalmanFittedStateOnPlane> kfsop (new genfit::KalmanFittedStateOnPlane(
			*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate())));
	// extrapolate back to reference plane.
	try {
		pathlenth = rep->extrapolateToPoint(*kfsop, P);
	} catch (genfit::Exception& e) {
		std::cerr << "Exception, next track" << std::endl;
		std::cerr << e.what();
		//delete kfsop;
		return WILD_DOUBLE;
	}

	state = *dynamic_cast<genfit::MeasuredStateOnPlane*> (kfsop.get());

	//delete kfsop;

	return pathlenth;
}

genfit::MeasuredStateOnPlane*  Track::extrapolateToPoint(TVector3 P, const int tr_point_id) const
{
	genfit::MeasuredStateOnPlane* state = new genfit::MeasuredStateOnPlane();
	double pathlenth = this->extrapolateToPoint(*state, P, tr_point_id);
	if(pathlenth <= WILD_DOUBLE) {
		delete state;
		return NULL;
	}
	else
		return state;
}

double Track::get_charge() const {
	double charge =  WILD_DOUBLE;

	genfit::AbsTrackRep* rep = _track->getCardinalRep();
	if(rep) {
		std::unique_ptr<genfit::StateOnPlane> state (this->extrapolateToLine(TVector3(0, 0, 0),
				TVector3(1, 0, 0)));
		if (state)
			charge = rep->getCharge(*state);
		//delete state;
	}

	return charge;
}

} //End of PHGenFit namespace
