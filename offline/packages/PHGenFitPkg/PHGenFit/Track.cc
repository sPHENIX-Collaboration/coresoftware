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

#define LogDebug(exp)		std::cout<<"DEBUG: "<<__FILE__<<": "<<__LINE__<<": "<< #exp <<" : "<< exp <<"\n"
#define LogError(exp)		std::cout<<"ERROR: "<<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"
#define LogWarning(exp)	std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"

#define WILD_DOUBLE -999999


namespace PHGenFit {

Track::Track(genfit::AbsTrackRep *rep, TVector3 seed_pos, TVector3 seed_mom, TMatrixDSym seed_cov)
{
//TODO Add input param check

	genfit::MeasuredStateOnPlane seedMSoP(rep);
	seedMSoP.setPosMomCov(seed_pos, seed_mom, seed_cov);
	const genfit::StateOnPlane seedSoP(seedMSoP);

	TVectorD seedState(6);
	TMatrixDSym seedCov(6);
	seedMSoP.get6DStateCov(seedState, seedCov);


	_track = new genfit::Track(rep, seedState, seedCov);
	//_track = NEW(genfit::Track)(rep, seedState, seedCov);
}

int Track::addMeasurements(std::vector<PHGenFit::Measurement*> &measurements)
{
	for(PHGenFit::Measurement* measurement : measurements)
	{
		std::vector<genfit::AbsMeasurement*> msmts;
		msmts.push_back(measurement->getMeasurement());
		_track->insertPoint(
				new genfit::TrackPoint(msmts, _track));

		_measurements.push_back(measurement);
	}

	//measurements.clear();

	return 0;
}

Track::~Track()
{
	delete _track;

	for(PHGenFit::Measurement* measurement : _measurements)
	{
		delete measurement;
	}
	_measurements.clear();
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

double Track::extrapolateToCylinder(genfit::MeasuredStateOnPlane& state, double radius, TVector3 line_point, TVector3 line_direction, const int tr_point_id) const
{
	double pathlenth = WILD_DOUBLE;

	genfit::AbsTrackRep* rep = _track->getCardinalRep();

//	genfit::TrackPoint* tp = _track->getPointWithMeasurementAndFitterInfo(
//			tr_point_id, rep);
//	if (tp == NULL) {
//		std::cout << "Track has no TrackPoint with fitterInfo! \n";
//		return WILD_DOUBLE;
//	}
	genfit::TrackPoint* tp = _track->getPointWithMeasurement(tr_point_id);
	if (tp == NULL) {
		std::cout << "tp == NULL! \n";
		return WILD_DOUBLE;
	}

	std::unique_ptr<genfit::MeasuredStateOnPlane> kfsop = NULL;
	if (dynamic_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))) {
		kfsop =
				std::unique_ptr < genfit::MeasuredStateOnPlane
						> (new genfit::KalmanFittedStateOnPlane(
								*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(
										rep))->getBackwardUpdate())));
	} else {
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
		std::cerr << "Exception, next track" << std::endl;
		std::cerr << e.what();
		//delete kfsop;
		return WILD_DOUBLE;
	}

	state = *dynamic_cast<genfit::MeasuredStateOnPlane*> (kfsop.get());

	//delete kfsop;

	return pathlenth;
}

genfit::MeasuredStateOnPlane*  Track::extrapolateToCylinder(double radius, TVector3 line_point, TVector3 line_direction, const int tr_point_id) const
{
	genfit::MeasuredStateOnPlane* state = new genfit::MeasuredStateOnPlane();
	double pathlenth = this->extrapolateToCylinder(*state, radius, line_point, line_direction);
	if(pathlenth <= WILD_DOUBLE) {
		delete state;
		return NULL;
	}
	else
		return state;
}

int Track::updateOneMeasurementKalman(
		const std::vector<PHGenFit::Measurement*>& measurements,
		std::map<double, PHGenFit::Track*>& incr_chi2s_new_tracks) const {

	//const int direction = 1;

	for (PHGenFit::Measurement* measurement : measurements) {

		PHGenFit::Track* new_track = new PHGenFit::Track(*this);

		genfit::Track *track = new_track->getGenFitTrack();
		genfit::AbsTrackRep* rep = track->getCardinalRep();
		genfit::TrackPoint *tp_base = track->getPointWithMeasurement(-1);

		//TODO use smart ptr
		std::unique_ptr<genfit::MeasuredStateOnPlane> state = NULL;
		genfit::SharedPlanePtr plane = NULL;

		bool newFi(!tp_base->hasFitterInfo(rep));
		if (newFi) {
			state = std::unique_ptr < genfit::MeasuredStateOnPlane
					> (new genfit::MeasuredStateOnPlane(rep));
			rep->setPosMomCov(*state, track->getStateSeed(),
					track->getCovSeed());
		} else {
			state =
					std::unique_ptr < genfit::MeasuredStateOnPlane
							> (new genfit::MeasuredStateOnPlane(
									static_cast<genfit::KalmanFitterInfo*>(tp_base->getFitterInfo(
											rep))->getFittedState(true)));
		}

		std::vector<genfit::AbsMeasurement*> msmts;
		msmts.push_back(measurement->getMeasurement());
		genfit::TrackPoint *tp = new genfit::TrackPoint(msmts, track);

		genfit::KalmanFitterInfo* fi = new genfit::KalmanFitterInfo(tp, rep);
		tp->setFitterInfo(fi);

		const std::vector<genfit::AbsMeasurement*>& rawMeasurements =
				tp->getRawMeasurements();
		// construct plane with first measurement
		plane = rawMeasurements[0]->constructPlane(*state);

		//double extLen = rep->extrapolateToPlane(*state, plane);
		rep->extrapolateToPlane(*state, plane);

		TVectorD stateVector(state->getState());
		TMatrixDSym cov(state->getCov());

		for (std::vector<genfit::AbsMeasurement*>::const_iterator it =
				rawMeasurements.begin(); it != rawMeasurements.end(); ++it) {
			fi->addMeasurementsOnPlane(
					(*it)->constructMeasurementsOnPlane(*state));
		}

		double chi2inc = 0;
		double ndfInc = 0;

		// update(s)
		const std::vector<genfit::MeasurementOnPlane *>& measurements_on_plane =
				fi->getMeasurementsOnPlane();
		for (std::vector<genfit::MeasurementOnPlane *>::const_iterator it =
				measurements_on_plane.begin();
				it != measurements_on_plane.end(); ++it) {
			const genfit::MeasurementOnPlane& mOnPlane = **it;
			//const double weight = mOnPlane.getWeight();

			const TVectorD& measurement(mOnPlane.getState());
			const genfit::AbsHMatrix* H(mOnPlane.getHMatrix());
			// (weighted) cov
			const TMatrixDSym& V(mOnPlane.getCov());

			TVectorD res(measurement - H->Hv(stateVector));

			// If hit, do Kalman algebra.
			{
				// calculate kalman gain ------------------------------
				// calculate covsum (V + HCH^T) and invert
				TMatrixDSym covSumInv(cov);
				H->HMHt(covSumInv);
				covSumInv += V;
				genfit::tools::invertMatrix(covSumInv);

				TMatrixD CHt(H->MHt(cov));
				TVectorD update(
						TMatrixD(CHt, TMatrixD::kMult, covSumInv) * res);
				//TMatrixD(CHt, TMatrixD::kMult, covSumInv).Print();

				stateVector += update;
				covSumInv.Similarity(CHt); // with (C H^T)^T = H C^T = H C  (C is symmetric)
				cov -= covSumInv;
			}

			TVectorD resNew(measurement - H->Hv(stateVector));

			// Calculate chi²
			TMatrixDSym HCHt(cov);
			H->HMHt(HCHt);
			HCHt -= V;
			HCHt *= -1;

			genfit::tools::invertMatrix(HCHt);

			chi2inc += HCHt.Similarity(resNew);

			ndfInc += measurement.GetNrows();
		} //loop measurements_on_plane

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
