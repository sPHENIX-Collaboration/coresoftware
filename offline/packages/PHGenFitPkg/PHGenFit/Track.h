/*!
 *  \file		Track.h
 *  \brief		Data structure and output of the fitting.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __PHGenFit_Track__
#define __PHGenFit_Track__

//STL
#include <vector>

//BOOST
//#include<boost/make_shared.hpp>
//
//#define SMART(expr) boost::shared_ptr<expr>
//#define NEW(expr) boost::make_shared<expr>

//GenFit


namespace genfit {

class AbsTrackRep;
class StateOnPlane;
class Track;

}

namespace PHGenFit {

class Measurement;

class Track
{
public:

	//! Default ctor
	Track(genfit::AbsTrackRep *rep, TVector3 seed_pos, TVector3 seed_mom, TMatrixDSym seed_cov);

	//! Default dtor
	~Track();

	//! Add measurement
	int addMeasurements(std::vector<PHGenFit::Measurement*> &measurements);

	/*!
	 * track_point 0 is the first one, and -1 is the last one
	 */
	double extrapolateToPlane(genfit::MeasuredStateOnPlane& state, TVector3 O, TVector3 n, const int tr_point_id = 0) const;
	//!
	double extrapolateToLine(genfit::MeasuredStateOnPlane& state, TVector3 line_point, TVector3 line_direction, const int tr_point_id = 0) const;
	//!
	double extrapolateToCylinder(genfit::MeasuredStateOnPlane& state, double radius, TVector3 line_point, TVector3 line_direction, const int tr_point_id = 0) const;
	//!
	double extrapolateToPoint(genfit::MeasuredStateOnPlane& state, TVector3 P, const int tr_point_id = 0) const;

	//!
	genfit::MeasuredStateOnPlane* extrapolateToPlane(TVector3 O, TVector3 n, const int tr_point_id = 0) const;
	//!
	genfit::MeasuredStateOnPlane* extrapolateToLine(TVector3 line_point, TVector3 line_direction, const int tr_point_id = 0) const;
	//!
	genfit::MeasuredStateOnPlane* extrapolateToCylinder(double radius, TVector3 line_point, TVector3 line_direction, const int tr_point_id = 0) const;
	//!
	genfit::MeasuredStateOnPlane* extrapolateToPoint(TVector3 P, const int tr_point_id = 0) const;
	//!
	genfit::Track* getGenFitTrack() {return _track;}

	double get_chi2() const {
		genfit::AbsTrackRep* rep = _track->getCardinalRep();
		double chi2 = _track->getFitStatus(rep)->getChi2();
		return chi2;
	}

	double get_ndf() const {
		genfit::AbsTrackRep* rep = _track->getCardinalRep();
		double ndf = _track->getFitStatus(rep)->getNdf();
		return ndf;
	}

	double get_charge() const;

	//SMART(genfit::Track) getGenFitTrack() {return _track;}

private:

	genfit::Track* _track;
	std::vector<PHGenFit::Measurement*> _measurements;
	//SMART(genfit::Track) _track;
};
} //End of PHGenFit namespace

#endif //__PHGenFit_Track__
