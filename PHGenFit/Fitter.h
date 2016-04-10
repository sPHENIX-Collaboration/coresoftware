#ifndef __PHGenFit_Fitter__
#define __PHGenFit_Fitter__

//STL
#include <vector>

//BOOST
//#include <boost/shared_ptr.hpp>
//#include <boost/make_shared.hpp>

//ROOT
#include <TDatabasePDG.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TRandom.h>
#include <TVector3.h>
#include "TDatabasePDG.h"
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>

//GenFit
#include <GenFit/ConstField.h>
#include <GenFit/Exception.h>
#include <GenFit/FieldManager.h>
#include <GenFit/KalmanFitterRefTrack.h>
#include <GenFit/StateOnPlane.h>
#include <GenFit/Track.h>
#include <GenFit/TrackPoint.h>

#include <GenFit/MaterialEffects.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/TGeoMaterialInterface.h>

#include <GenFit/EventDisplay.h>

#include <GenFit/HelixTrackModel.h>
#include <GenFit/MeasurementCreator.h>

#include "GenFit/PlanarMeasurement.h"
#include "GenFit/DetPlane.h"
#include "GenFit/SharedPlanePtr.h"

#include <GenFit/KalmanFittedStateOnPlane.h>
//#include <AbsKalmanFitter.h>
//#include <KalmanFitter.h>
//#include <KalmanFitterRefTrack.h>
#include <GenFit/KalmanFitterInfo.h>
//#include <KalmanFitStatus.h>
#include <GenFit/KalmanFitter.h>

//GenFitExp
#include <genfitexp/Field2D.h>

//PHGenFit
#include <Track.h>

#define LogDEBUG    std::cout<<"DEBUG: "<<__LINE__<<"\n"

namespace PHGenFit {

class Fitter
{
public:
	//! Default constructor
	Fitter(const std::string tgeo_file_name,
			const std::string field_file_name,
			const double field_scaling_factor = 1.4/1.5,
			const std::string fitter_choice = "KalmanFitterRefTrack",
			const std::string track_rep_choice = "RKTrackRep");
	//! Default destructor
	~Fitter();

	int processTrack(PHGenFit::Track* track);

private:

	TGeoManager* _tgeo_manager;

	genfit::EventDisplay* _display;
	genfit::AbsKalmanFitter* _fitter;







}; //class Fitter

} //End of PHGenFit namespace

#endif //__PHGenFit_Fitter__
