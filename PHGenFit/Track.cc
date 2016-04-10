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


namespace PHGenFit {

Track::Track(genfit::AbsTrackRep *rep, TVector3 seed_pos, TVector3 seed_mom, TMatrixDSym seed_cov)
{
//FIXME Add input param check

	genfit::MeasuredStateOnPlane seedMSoP(rep);
	seedMSoP.setPosMomCov(seed_pos, seed_mom, seed_cov);
	const genfit::StateOnPlane seedSoP(seedMSoP);

	TVectorD seedState(6);
	TMatrixDSym seedCov(6);
	seedMSoP.get6DStateCov(seedState, seedCov);


	_track = new genfit::Track(rep, seedState, seedCov);
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

genfit::StateOnPlane* Track::extrapolateToLine(TVector3 line_point, TVector3 line_direction) const
{
	genfit::AbsTrackRep* rep = _track->getCardinalRep();
	genfit::TrackPoint* tp = _track->getPointWithMeasurementAndFitterInfo(
			0, rep);
	if (tp == NULL) {
		std::cout << "Track has no TrackPoint with fitterInfo! \n";
		return NULL;
	}
	genfit::KalmanFittedStateOnPlane *kfsop = new genfit::KalmanFittedStateOnPlane(
			*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()));
	// extrapolate back to reference plane.
	try {
		rep->extrapolateToLine(*kfsop, line_point, line_direction);
	} catch (genfit::Exception& e) {
		std::cerr << "Exception, next track" << std::endl;
		std::cerr << e.what();
		return NULL;
	}

	return kfsop;
}

} //End of PHGenFit namespace
