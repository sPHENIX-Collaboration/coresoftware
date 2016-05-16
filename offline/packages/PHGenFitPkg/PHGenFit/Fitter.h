/*!
 *  \file		Fitter.h
 *  \brief		Fitter class handles setups for the fitting.
 *  \details	Fitter class handles setups for the fitting like Geometry, Fields, fitter choice, etc.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __PHGenFit_Fitter__
#define __PHGenFit_Fitter__

//STL
#include <vector>

//BOOST
//#include<boost/make_shared.hpp>
//
//#define SMART(expr) boost::shared_ptr<expr>
//#define NEW(expr) boost::make_shared<expr>

class TGeoManager;

namespace genfit{
	class FieldManager;
	class MaterialEffects;
	class EventDisplay;
	class AbsKalmanFitter;
}

namespace PHGenFit {
class Track;

class Fitter
{
public:
	//! Default constructor
	Fitter(const std::string tgeo_file_name,
			const std::string field_file_name,
			const double field_scaling_factor = 1.4/1.5,
			const std::string fitter_choice = "KalmanFitterRefTrack",
			const std::string track_rep_choice = "RKTrackRep",
			const bool doEventDisplay = false);
	//! Default destructor
	~Fitter();

	int processTrack(PHGenFit::Track* track, const bool save_to_evt_disp = false);

	int displayEvent();

private:

	bool _doEventDisplay;

	TGeoManager* _tgeo_manager;

	genfit::EventDisplay* _display;
	genfit::AbsKalmanFitter* _fitter;







}; //class Fitter

} //End of PHGenFit namespace

#endif //__PHGenFit_Fitter__
