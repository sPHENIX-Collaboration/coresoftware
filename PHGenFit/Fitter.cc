#include <Fitter.h>

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
	delete _tgeo_manager;
	delete _display;
	delete _fitter;
}

int Fitter::processTrack(PHGenFit::Track* track)
{

//FIXME Add more fitting info
	_fitter->processTrack(track->getGenFitTrack());

	return 0;
}

} //End of PHGenFit namespace
