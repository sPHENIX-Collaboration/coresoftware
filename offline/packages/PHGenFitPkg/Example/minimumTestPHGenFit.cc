/*!
 *  \file		minimumTestPHGenFit.cc
 *  \brief		Minimum program to demonstrate the usage of PHGenFit.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

//STL
#include <vector>

//ROOT
#include <TVector3.h>
#include <TMatrixDSym.h>

//GenFit
#include <GenFit/AbsTrackRep.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/StateOnPlane.h>

//PHGenFit
#include <phgenfit/Fitter.h>
#include <phgenfit/Track.h>
#include <phgenfit/Measurement.h>
#include <phgenfit/PlanarMeasurement.h>
#include <phgenfit/SpacepointMeasurement.h>

#include <phfield/PHFieldUtility.h>

#define LogDEBUG    std::cout<<"DEBUG: "<<__LINE__<<"\n"

void get_seed(TVector3& seed_pos, TVector3& seed_mom, TMatrixDSym& seed_cov)
{
	seed_pos.SetXYZ(0,0,0);
	seed_mom.SetXYZ(10,-5,0);
	seed_cov.ResizeTo(6,6);
}

std::vector<TVector3> get_raw_measurements()
{
	std::vector<TVector3> v_pos;
	v_pos.push_back(TVector3(2.22459,-1.54767,-2.37792));
	v_pos.push_back(TVector3(3.80050,-2.64444,-2.16561));
	v_pos.push_back(TVector3(7.80344,-5.41815,-1.98928));
	v_pos.push_back(TVector3(8.63214,-5.97797,-1.59626));
	return v_pos;
}

int main(int argc, char**argv) {

	//! Initiallize Geometry, Field, Fitter
	PHGenFit::Fitter* fitter = new PHGenFit::Fitter("sPHENIX_Geo.root",
	    PHFieldUtility::BuildFieldMap(PHFieldUtility::DefaultFieldConfig(), 1),
	    "KalmanFitter","RKTrackRep",false);

	//! Build TrackRep from particle assumption
	int pid = -13; //mu+
	genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pid);

	//! Initiallize track with seed from pattern recognition
	TVector3 seed_pos;
	TVector3 seed_mom;
	TMatrixDSym seed_cov;
	get_seed(seed_pos,seed_mom, seed_cov);
	PHGenFit::Track* track = new PHGenFit::Track(rep, seed_pos,seed_mom, seed_cov);

	//! Create measurements
	std::vector<TVector3> v_pos = get_raw_measurements();
	double res_phi = 0.005; //cm
	double res_z = 0.04; //cm
	std::vector<PHGenFit::Measurement*> measurements;
	for (unsigned int imeasurement = 0; imeasurement < v_pos.size(); imeasurement++) {
		TVector3 pos = v_pos[imeasurement];
		TVector3 n(pos.x(),pos.Y(),0);
		PHGenFit::Measurement* meas = new PHGenFit::PlanarMeasurement(pos,n,res_phi, res_z);
		//PHGenFit::Measurement* meas = new PHGenFit::SpacepointMeasurement(pos,res_phi);
		meas->getMeasurement()->Print();
		measurements.push_back(meas);
	}

	//! Add measurements to track
	track->addMeasurements(measurements);

	//! Fit the track
	fitter->processTrack(track, false);

	//! Extrapolate to beam line
	//genfit::MeasuredStateOnPlane* state = track->extrapolateToLine(TVector3(0, 0, 0), TVector3(0, 0, 1));
	//genfit::MeasuredStateOnPlane* state = track->extrapolateToPoint(TVector3(0, 0, 0));
	genfit::MeasuredStateOnPlane* state = track->extrapolateToCylinder(1.,TVector3(0, 0, 0), TVector3(0, 0, 1));
	state->Print();
	delete state;

	//! Event display, uncomment to use
	//fitter->displayEvent();

	//! Comment off if want to keep event display.

	delete track;

	delete fitter;

	return 0;
}
