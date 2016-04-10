/*!
 *  \file		testPHGenFit.cc
 *  \brief		Program to demonstrate the usage of PHGenFit.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include <vector>

#include <TVector3.h>
#include <TMatrixDSym.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TRandom.h>

#include <GenFit/AbsTrackRep.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/StateOnPlane.h>

#include <phgenfit/Fitter.h>
#include <phgenfit/Track.h>
#include <phgenfit/Measurement.h>
#include <phgenfit/PlanarMeasurement.h>

#define LogDEBUG    std::cout<<"DEBUG: "<<__LINE__<<"\n"

int main(int argc, char**argv) {
	//! Initiallize Geometry, Field, Fitter
	PHGenFit::Fitter* fitter = new PHGenFit::Fitter("sPHENIX_Geo.root",
			"sPHENIX.2d.root", 1.4 / 1.5);

	TFile *fPHG4Hits = TFile::Open("AnaSvtxTracksForGenFit.root", "read");
	if (!fPHG4Hits) {
		std::cout << "No TFile Openned: " << __LINE__ << "\n";
		return -1;
	}
	TTree *T = (TTree*) fPHG4Hits->Get("tracks");
	if (!T) {
		std::cout << "No TTree Found: " << __LINE__ << "\n";
		return -1;
	}

	double resolution_detector_xy = 0.005/3.; //50/3. micron

#define NLAYERS 7


	Float_t Cluster_x[NLAYERS];
	Float_t Cluster_y[NLAYERS];
	Float_t Cluster_z[NLAYERS];
	Float_t Cluster_size_dphi[NLAYERS];
	Float_t Cluster_size_dz[NLAYERS];
	Float_t True_px;
	Float_t True_py;
	Float_t True_pz;
	Float_t True_vx;
	Float_t True_vy;
	Float_t True_vz;
	Float_t AlanDion_px;
	Float_t AlanDion_py;
	Float_t AlanDion_pz;
	Float_t AlanDion_dca2d;
	Int_t nhits;

	T->SetBranchAddress("nhits", &nhits);
	T->SetBranchAddress("gpx", &True_px);
	T->SetBranchAddress("gpy", &True_py);
	T->SetBranchAddress("gpz", &True_pz);
	T->SetBranchAddress("gvx", &True_vx);
	T->SetBranchAddress("gvy", &True_vy);
	T->SetBranchAddress("gvz", &True_vz);
	T->SetBranchAddress("px", &AlanDion_px);
	T->SetBranchAddress("py", &AlanDion_py);
	T->SetBranchAddress("pz", &AlanDion_pz);
	T->SetBranchAddress("dca2d", &AlanDion_dca2d);
	T->SetBranchAddress("x", Cluster_x);
	T->SetBranchAddress("y", Cluster_y);
	T->SetBranchAddress("z", Cluster_z);
	T->SetBranchAddress("size_dphi", Cluster_size_dphi);
	T->SetBranchAddress("size_dz", Cluster_size_dz);

	for (unsigned int ientry = 0; ientry < 1; ++ientry) {
		//T->GetEntry(atoi(argv[1]));
		T->GetEntry(ientry);

		if (nhits < 0) {
			LogDEBUG;
			continue;
		}

		// true start values
		TVector3 init_pos(0, 0, 0); //cm
		TVector3 True_mom(True_px, True_py, True_pz);

		// Seed: use smeared values
		const bool smearPosMom = true; // init the Reps with smeared init_pos and True_mom
		const double posSmear = 10 * resolution_detector_xy;     // cm
		const double momSmear = 3. / 180. * TMath::Pi();     // rad
		const double momMagSmear = 0.1;   // relative

		TVector3 seed_pos(init_pos);
		TVector3 seed_mom(True_mom);
		if (smearPosMom) {
			seed_pos.SetX(gRandom->Gaus(seed_pos.X(), posSmear));
			seed_pos.SetY(gRandom->Gaus(seed_pos.Y(), posSmear));
			seed_pos.SetZ(gRandom->Gaus(seed_pos.Z(), posSmear));

			seed_mom.SetPhi(gRandom->Gaus(True_mom.Phi(), momSmear));
			seed_mom.SetTheta(gRandom->Gaus(True_mom.Theta(), momSmear));
			seed_mom.SetMag(
					gRandom->Gaus(True_mom.Mag(),
							momMagSmear * True_mom.Mag()));
		}

		// approximate covariance
		TMatrixDSym seed_cov(6);

		for (int idim = 0; idim < 3; ++idim)
			seed_cov(idim, idim) = resolution_detector_xy
					* resolution_detector_xy;
		for (int idim = 3; idim < 6; ++idim)
			seed_cov(idim, idim) = pow(
					resolution_detector_xy / NLAYERS / sqrt(3), 2);


		//! Build TrackRep from particle assumption
		int pid = -13; //mu+
		genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pid);

		//! Initiallize track with seed from pattern recognition
		PHGenFit::Track* track = new PHGenFit::Track(rep, seed_pos,
				seed_mom, seed_cov);

		//! Create measurements
		std::vector<PHGenFit::Measurement*> measurements;
		for (int imeas = 0; imeas < NLAYERS; imeas++) {
			TVector3 pos(Cluster_x[imeas],Cluster_y[imeas],Cluster_z[imeas]);
			TVector3 n(Cluster_x[imeas],Cluster_y[imeas],0);
			TVector3 v(0,0,1);
			TVector3 u = v.Cross(n);

			PHGenFit::Measurement* meas = new PHGenFit::PlanarMeasurement(
					pos,u,v,
					Cluster_size_dphi[imeas],Cluster_size_dz[imeas]);
			measurements.push_back(meas);
		}

		//! Add measurements to track
		track->addMeasurements(measurements);

		//! Fit the track
		fitter->processTrack(track, true);

		//!
		genfit::StateOnPlane* state_at_beam_line = track->extrapolateToLine(
				TVector3(0, 0, 0), TVector3(0, 0, 1));
		state_at_beam_line->Print();
	}


	//! Event display
	fitter->displayEvent();

	return 0;
}
