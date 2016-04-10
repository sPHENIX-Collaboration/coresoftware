#ifndef __PHGenFit_Track__
#define __PHGenFit_Track__

#include <vector>

#include <GenFit/Track.h>


#include <Measurement.h>

namespace PHGenFit {

class Track
{
public:

	//! Default ctor
	Track(genfit::AbsTrackRep *rep, TVector3 seed_pos, TVector3 seed_mom, TMatrixDSym seed_cov);

	//! Default dtor
	~Track();

	int addMeasurements(std::vector<PHGenFit::Measurement*> measurements);

	//!
	genfit::Track* getGenFitTrack() {return _track;}

private:
	genfit::Track* _track;
};
} //End of PHGenFit namespace

#endif //__PHGenFit_Track__
