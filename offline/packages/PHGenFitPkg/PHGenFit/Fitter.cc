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
#include <GenFit/DAF.h>
#include <GenFit/RKTrackRep.h>

//GenFitExp
#include <genfitexp/Field2D.h>

//PHGenFit
#include "Fitter.h"
#include "Track.h"

#define LogDEBUG(exp)		std::cout<<"DEBUG: "<<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"
#define LogERROR(exp)		std::cout<<"ERROR: "<<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"
#define LogWARNING(exp)	std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp <<"\n"

namespace PHGenFit {

Fitter::Fitter(
		const std::string tgeo_file_name,
		const std::string field_file_name,
		const double field_scaling_factor,
		const std::string fitter_choice,
		const std::string track_rep_choice,
		const bool doEventDisplay
) : verbosity(0), _doEventDisplay(doEventDisplay)
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
	else if(fitter_choice.compare("DafSimple")==0)
		_fitter = new genfit::DAF(false);
	else if(fitter_choice.compare("DafRef")==0)
		_fitter = new genfit::DAF(true);
	else
		_fitter = new genfit::KalmanFitter();

	genfit::Exception::quiet(true);
}

Fitter::~Fitter()
{
	if(_fitter)
		delete _fitter;
	if(_tgeo_manager)
		//delete _tgeo_manager;
		_tgeo_manager->Delete();
	if(_display)
		delete _display;
}

int Fitter::processTrack(PHGenFit::Track* track, const bool save_to_evt_disp) {
	genfit::Track* fitTrack = track->getGenFitTrack();
	if(!fitTrack->checkConsistency()){
		if(verbosity >= 2) LogWARNING("genfit::Track::checkConsistency() failed!") ;
		return -1;
	}

	try {
		_fitter->processTrack(fitTrack);
	} catch (genfit::Exception& e) {
		if (verbosity >= 1) {
			std::cerr << "PHGenFit::Exception: \n";
			std::cerr << e.what();
			std::cerr << "Exception, next track" << std::endl;
		}
		return -1;
	}

	if(!fitTrack->checkConsistency()){
		if(verbosity >= 2) LogWARNING("genfit::Track::checkConsistency() failed!") ;
		return -1;
	}

	genfit::AbsTrackRep* rep = fitTrack->getCardinalRep();
	if (!fitTrack->getFitStatus(rep)->isFitConverged()) {
		if(verbosity >= 2) LogWARNING("Track could not be fitted successfully! Fit is not converged!");
		return -1;
	}

	if (_display and save_to_evt_disp)
		_display->addEvent(track->getGenFitTrack());

	return 0;
}

Fitter* Fitter::getInstance(const std::string tgeo_file_name,
		const std::string field_file_name, const double field_scaling_factor,
		const std::string fitter_choice, const std::string track_rep_choice,
		const bool doEventDisplay) {

	TGeoManager* tgeo_manager = TGeoManager::Import(tgeo_file_name.data(), "Default");
	if(!tgeo_manager)
	{
		LogERROR("No TGeoManager found!");
		return NULL;
	}

	genfit::Field2D *fieldMap = new genfit::Field2D();
	if(!fieldMap->initialize(field_file_name.data()))
	{
		LogERROR("Field map initialization failed!");
		delete fieldMap;
		return NULL;
	}
	fieldMap->re_scale(field_scaling_factor);// Re-scale to 1.4 T

	return new Fitter(tgeo_manager, fieldMap, fitter_choice, track_rep_choice, doEventDisplay);
}

Fitter::Fitter(TGeoManager* tgeo_manager, genfit::AbsBField* fieldMap,
		const std::string fitter_choice, const std::string track_rep_choice,
		const bool doEventDisplay): verbosity(0), _tgeo_manager(tgeo_manager), _doEventDisplay(doEventDisplay)
{

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

int Fitter::displayEvent()
{
	if(_display)
		_display->open();
	else
		if(verbosity >= 0) LogERROR("No genfit::EventDisplay found!");

	return 0;
}
} //End of PHGenFit namespace
