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
#include<boost/make_shared.hpp>

#define SMART(expr) boost::shared_ptr<expr>
#define NEW(expr) boost::make_shared<expr>

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
	int addMeasurements(std::vector<PHGenFit::Measurement*> measurements);

	//!
	genfit::StateOnPlane* extrapolateToLine(TVector3 line_point, TVector3 line_direction) const;

	//!
	genfit::Track* getGenFitTrack() {return _track.get();}
	//SMART(genfit::Track) getGenFitTrack() {return _track;}

private:

	//genfit::Track* _track;
	SMART(genfit::Track) _track;

//	TODO figure out how to handle multiple TrackReps
//	double _chi2;
//	double _ndf;

};
} //End of PHGenFit namespace

#endif //__PHGenFit_Track__
