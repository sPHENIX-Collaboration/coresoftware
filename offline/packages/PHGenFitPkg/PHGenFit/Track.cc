/*!
 *  \file		Track.cc
 *  \brief		Data structure and output of the fitting.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

//STL

//BOOST
#include <boost/foreach.hpp>

//GenFit
#include <GenFit/Track.h>
#include <GenFit/KalmanFittedStateOnPlane.h>
#include <GenFit/KalmanFitterInfo.h>
#include <GenFit/KalmanFitter.h>

//PHGenFit
#include "Track.h"
#include "Measurement.h"

#define LogDebug(exp)		std::cout<<"DEBUG: "<<__FILE__<<": "<<__LINE__<<": "<< #exp <<" : "<< exp <<"\n"
#define LogError(exp)		std::cout<<"ERROR: "<<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"
#define LogWarning(exp)	std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"


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

int Track::addMeasurements(std::vector<PHGenFit::Measurement*> measurements)
{
	BOOST_FOREACH(PHGenFit::Measurement* measurement, measurements)
	{
		std::vector<genfit::AbsMeasurement*> msmts;
		msmts.push_back(measurement->getMeasurement());
		_track->insertPoint(
				new genfit::TrackPoint(msmts, _track));
	}

	return 0;
}

Track::~Track()
{
	delete _track;
}

double Track::extrapolateToPlane(genfit::MeasuredStateOnPlane* state, TVector3 O, TVector3 n, const int tr_point_id) const
{
	double pathlenth = -9999;

	genfit::SharedPlanePtr destPlane(new genfit::DetPlane(O, n));

	genfit::AbsTrackRep* rep = _track->getCardinalRep();
	genfit::TrackPoint* tp = _track->getPointWithMeasurementAndFitterInfo(
			tr_point_id, rep);
	if (tp == NULL) {
		std::cout << "Track has no TrackPoint with fitterInfo! \n";
		return pathlenth;
	}

	state = new genfit::KalmanFittedStateOnPlane(
			*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()));
	// extrapolate back to reference plane.
	try {
		pathlenth = rep->extrapolateToPlane(*state, destPlane);
	} catch (genfit::Exception& e) {
		std::cerr << "Exception, next track" << std::endl;
		std::cerr << e.what();
		return -9999;
	}

	return pathlenth;
}

genfit::MeasuredStateOnPlane* Track::extrapolateToPlane(TVector3 O, TVector3 n, const int tr_point_id) const
{
	genfit::MeasuredStateOnPlane* state = NULL;
	double pathlenth = this->extrapolateToPlane(state, O, n, tr_point_id);
	if (pathlenth > 0)
		return state;
	else
		return NULL;
}

double Track::extrapolateToLine(genfit::MeasuredStateOnPlane* state, TVector3 line_point, TVector3 line_direction, const int tr_point_id) const
{
	double pathlenth = -9999;

	genfit::AbsTrackRep* rep = _track->getCardinalRep();
	genfit::TrackPoint* tp = _track->getPointWithMeasurementAndFitterInfo(
			tr_point_id, rep);
	if (tp == NULL) {
		std::cout << "Track has no TrackPoint with fitterInfo! \n";
		return -9999;
	}
	state = new genfit::KalmanFittedStateOnPlane(
			*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()));
	// extrapolate back to reference plane.
	try {
		pathlenth = rep->extrapolateToLine(*state, line_point, line_direction);
	} catch (genfit::Exception& e) {
		std::cerr << "Exception, next track" << std::endl;
		std::cerr << e.what();
		return -9999;
	}

	return pathlenth;
}

genfit::MeasuredStateOnPlane* Track::extrapolateToLine(TVector3 line_point, TVector3 line_direction, const int tr_point_id) const
{
	genfit::MeasuredStateOnPlane* state = NULL;
	double pathlenth = this->extrapolateToLine(state, line_point, line_direction, tr_point_id);
	if (pathlenth > 0)
		return state;
	else
		return NULL;
}

double Track::extrapolateToCylinder(genfit::MeasuredStateOnPlane* state, double radius, TVector3 line_point, TVector3 line_direction, const int tr_point_id) const
{
	double pathlenth = -9999;

	genfit::AbsTrackRep* rep = _track->getCardinalRep();
	genfit::TrackPoint* tp = _track->getPointWithMeasurementAndFitterInfo(
			tr_point_id, rep);
	if (tp == NULL) {
		std::cout << "Track has no TrackPoint with fitterInfo! \n";
		return -9999;
	}
	state = new genfit::KalmanFittedStateOnPlane(
			*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getForwardUpdate()));
	// extrapolate back to reference plane.
	try {
		//rep->extrapolateToLine(*kfsop, line_point, line_direction);
		pathlenth = rep->extrapolateToCylinder(*state, radius, line_point, line_direction);
	} catch (genfit::Exception& e) {
		std::cerr << "Exception, next track" << std::endl;
		std::cerr << e.what();
		return -9999;
	}

	return pathlenth;
}

genfit::MeasuredStateOnPlane*  Track::extrapolateToCylinder(double radius, TVector3 line_point, TVector3 line_direction, const int tr_point_id) const
{
	genfit::MeasuredStateOnPlane* state = NULL;
	double pathlenth = this->extrapolateToCylinder(state, radius, line_point, line_direction);
	if (pathlenth > 0)
		return state;
	else
		return NULL;
}

double Track::extrapolateToPoint(genfit::MeasuredStateOnPlane* state, TVector3 P, const int tr_point_id) const
{
	double pathlenth = -9999;
	genfit::AbsTrackRep* rep = _track->getCardinalRep();
	genfit::TrackPoint* tp = _track->getPointWithMeasurementAndFitterInfo(
			tr_point_id, rep);
	if (tp == NULL) {
		std::cout << "Track has no TrackPoint with fitterInfo! \n";
		return -9999;
	}
	state = new genfit::KalmanFittedStateOnPlane(
			*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()));
	// extrapolate back to reference plane.
	try {
		pathlenth = rep->extrapolateToPoint(*state, P);
	} catch (genfit::Exception& e) {
		std::cerr << "Exception, next track" << std::endl;
		std::cerr << e.what();
		return -9999;
	}

	return pathlenth;
}

genfit::MeasuredStateOnPlane*  Track::extrapolateToPoint(TVector3 P, const int tr_point_id) const
{
	genfit::MeasuredStateOnPlane* state = NULL;
	double pathlenth = this->extrapolateToPoint(state, P, tr_point_id);
	if (pathlenth > 0)
		return state;
	else
		return NULL;
}
} //End of PHGenFit namespace
