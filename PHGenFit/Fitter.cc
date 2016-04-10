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
//#include <TDatabasePDG.h>
//#include <TEveManager.h>
//#include <TGeoManager.h>
//#include <TRandom.h>
//#include <TVector3.h>
//#include "TDatabasePDG.h"
//#include <TMath.h>
//#include <TF1.h>
//#include <TH1D.h>
//#include <TH2D.h>
//#include <TProfile.h>
//#include <TFile.h>
//#include <TTree.h>
//#include <TCanvas.h>
//#include <TROOT.h>
//#include <TStyle.h>

//GenFit
//#include <GenFit/ConstField.h>
//#include <GenFit/Exception.h>
//#include <GenFit/FieldManager.h>
////#include <GenFit/KalmanFitterRefTrack.h>
////#include <GenFit/StateOnPlane.h>
////#include <GenFit/Track.h>
////#include <GenFit/TrackPoint.h>
////
//#include <GenFit/MaterialEffects.h>
//#include <GenFit/RKTrackRep.h>
//#include <GenFit/TGeoMaterialInterface.h>
////
//#include <GenFit/EventDisplay.h>
////
////#include <GenFit/HelixTrackModel.h>
////#include <GenFit/MeasurementCreator.h>
////
////#include "GenFit/PlanarMeasurement.h"
////#include "GenFit/DetPlane.h"
////#include "GenFit/SharedPlanePtr.h"
//
////#include <GenFit/KalmanFittedStateOnPlane.h>
//#include <AbsKalmanFitter.h>
//#include <KalmanFitter.h>
//#include <KalmanFitterRefTrack.h>
//#include <GenFit/KalmanFitterInfo.h>
//#include <GenFit/KalmanFitter.h>

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
		const std::string track_rep_choice
)
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
	_display = genfit::EventDisplay::getInstance();

	// init fitter
	_fitter = new genfit::KalmanFitterRefTrack();
}

Fitter::~Fitter()
{
	delete _fitter;
	delete _tgeo_manager;
	delete _display;
}

int Fitter::processTrack(PHGenFit::Track* track, const bool save_to_evt_disp)
{

//TODO Add savety checks
	_fitter->processTrack(track->getGenFitTrack());

	if(save_to_evt_disp)
		_display->addEvent(track->getGenFitTrack());

	return 0;
}

int Fitter::displayEvent()
{
	_display->open();

	return 0;
}
} //End of PHGenFit namespace
