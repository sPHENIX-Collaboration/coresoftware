#ifndef __PHGenFit_Fitter__
#define __PHGenFit_Fitter__

#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>

#include <MaterialEffects.h>
#include <RKTrackRep.h>
//#include <G4eTrackRep.h>
#include <TGeoMaterialInterface.h>

#include <EventDisplay.h>

#include <HelixTrackModel.h>
#include <MeasurementCreator.h>

#include <TDatabasePDG.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TRandom.h>
#include <TVector3.h>
#include <vector>

#include "TDatabasePDG.h"
#include <TMath.h>

#include "PlanarMeasurement.h"
#include "DetPlane.h"
#include "SharedPlanePtr.h"
//#include <boost/shared_ptr.hpp>
//#include <boost/make_shared.hpp>
#include <KalmanFittedStateOnPlane.h>
//#include <AbsKalmanFitter.h>
//#include <KalmanFitter.h>
//#include <KalmanFitterRefTrack.h>
#include <KalmanFitterInfo.h>
//#include <KalmanFitStatus.h>
#include <KalmanFitter.h>


#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>

#include <Field2D.h>

#define LogDEBUG    std::cout<<"DEBUG: "<<__LINE__<<"\n"

namespace PHGenFit {

class Fitter
{
	};
} //End of PHGenFit namespace

#endif //__PHGenFit_Fitter__
