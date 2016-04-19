/*!
 *  \file		Fitter.cc
 *  \brief		Fitter class handles setups for the fitting.
 *  \details	Fitter class handles setups for the fitting like Geometry, Fields, fitter choice, etc.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

//BOOST
//#include <boost/shared_ptr.hpp>
//#include <boost/make_shared.hpp>

//ROOT
#include <TGeoManager.h>

//GenFit
#include <GenFit/FieldManager.h>
#include <GenFit/MaterialEffects.h>
#include <GenFit/TGeoMaterialInterface.h>
#include <GenFit/EventDisplay.h>
#include <GenFit/AbsKalmanFitter.h>
#include <GenFit/KalmanFitter.h>
#include <GenFit/KalmanFitterRefTrack.h>
#include <GenFit/RKTrackRep.h>

//GenFitExp
#include <genfitexp/Field2D.h>

//PHGenFit
#include "Fitter.h"
#include "Track.h"

#define LogDEBUG    std::cout<<"DEBUG: "<<__LINE__<<"\n"

namespace PHGenFit {

Fitter::Fitter(
		const std::string tgeo_file_name,
		const std::string field_file_name,
		const double field_scaling_factor,
		const std::string fitter_choice,
		const std::string track_rep_choice,
		const bool doEventDisplay
) : _doEventDisplay(doEventDisplay)
{
	_tgeo_manager = new TGeoManager("Default", "Geane geometry");
	TGeoManager::Import(tgeo_file_name.data());

	genfit::Field2D *fieldMap = new genfit::Field2D(field_file_name.data());
	fieldMap->re_scale(field_scaling_factor);// Re-scale to 1.4 T

	genfit::FieldManager::getInstance()->init(
			fieldMap);
	genfit::MaterialEffects::getInstance()->init(
			new genfit::TGeoMaterialInterface());

	// init event display
	if(_doEventDisplay)
		_display = genfit::EventDisplay::getInstance();
	else
		_display = NULL;

	// init fitter
	if(fitter_choice.compare("KalmanFitterRefTrack")==0)
		_fitter = new genfit::KalmanFitterRefTrack();
	else if(fitter_choice.compare("KalmanFitter")==0)
		_fitter = new genfit::KalmanFitter();
	else
		_fitter = new genfit::KalmanFitter();
}

Fitter::~Fitter()
{
	if(_fitter)
		delete _fitter;
	if(_tgeo_manager)
		delete _tgeo_manager;
	if(_display)
		delete _display;
}

int Fitter::processTrack(PHGenFit::Track* track, const bool save_to_evt_disp)
{

//TODO Add safety checks
	_fitter->processTrack(track->getGenFitTrack());

	if(_display and save_to_evt_disp)
		_display->addEvent(track->getGenFitTrack());

	return 0;
}

int Fitter::displayEvent()
{
	_display->open();

	return 0;
}
} //End of PHGenFit namespace
