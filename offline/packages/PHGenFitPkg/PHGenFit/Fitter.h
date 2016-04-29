/*!
 *  \file		Fitter.h
 *  \brief		Fitter class handles setups for the fitting.
 *  \details	Fitter class handles setups for the fitting like Geometry, Fields, fitter choice, etc.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __PHGenFit_Fitter__
#define __PHGenFit_Fitter__

//STL

#include <GenFit/EventDisplay.h>
#include <string>

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
	class AbsBField;
	class Field2D;
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

	Fitter(TGeoManager* tgeo_manager,
			genfit::AbsBField* fieldMap,
			const std::string fitter_choice = "KalmanFitterRefTrack",
			const std::string track_rep_choice = "RKTrackRep",
			const bool doEventDisplay = false);

	//! Default destructor
	~Fitter();

	static Fitter* getInstance(const std::string tgeo_file_name,
			const std::string field_file_name,
			const double field_scaling_factor = 1.4/1.5,
			const std::string fitter_choice = "KalmanFitterRefTrack",
			const std::string track_rep_choice = "RKTrackRep",
			const bool doEventDisplay = false);

	int processTrack(PHGenFit::Track* track, const bool save_to_evt_disp = false);

	int displayEvent();

	bool is_do_Event_Display() const {
		return _doEventDisplay;
	}

	void set_do_Event_Display(bool doEventDisplay) {
		_doEventDisplay = doEventDisplay;
		if(!_display && _doEventDisplay)
			_display = genfit::EventDisplay::getInstance();
	}

	genfit::EventDisplay* getEventDisplay()
	{
		return _display;
	}

private:

	TGeoManager* _tgeo_manager;

	bool _doEventDisplay;

	genfit::EventDisplay* _display;
	genfit::AbsKalmanFitter* _fitter;







}; //class Fitter

} //End of PHGenFit namespace

#endif //__PHGenFit_Fitter__
