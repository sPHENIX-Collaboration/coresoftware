/*!
 *  \file		testPHGenFit.cc
 *  \brief		Program to demonstrate the usage of PHGenFit.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

//STL
#include <vector>

//BOOST
#include<boost/make_shared.hpp>

#define SMART(expr) boost::shared_ptr<expr>
#define NEW(expr) boost::make_shared<expr>

#include <phfield/PHFieldUtility.h>

//ROOT
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

//GenFit
#include <GenFit/AbsTrackRep.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/StateOnPlane.h>

//PHGenFit
#include <phgenfit/Fitter.h>
#include <phgenfit/Track.h>
#include <phgenfit/Measurement.h>
#include <phgenfit/PlanarMeasurement.h>

#define LogDEBUG    std::cout<<"DEBUG: "<<__LINE__<<"\n"

//void pause() {
//  std::cout << "Press ENTER to continue..." << std::flush;
//  std::cin.clear();  // use only if a human is involved
//  std::cin.flush();  // use only if a human is involved
//  std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
//}

int main(int argc, char**argv) {

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

	//! Initiallize Geometry, Field, Fitter
	PHGenFit::Fitter* fitter = new PHGenFit::Fitter("sPHENIX_Geo.root",PHFieldUtility::BuildFieldMap(PHFieldUtility::DefaultFieldConfig(), 1));


	double resolution_detector_xy = 0.005/3.; //50/3. micron


	TH2D *hpT_residual_vs_pT = new TH2D("hpT_residual_vs_pT", "#Delta pT/pT; pT[GeV/c]; #Delta pT/pT", 40, 0.5, 40.5, 1000, -1, 1);

	TH2D *hDCAr_vs_pT = new TH2D("hDCAr_vs_pT", "DCAr vs. p; p [GeV/c]; DCAr [cm]", 40, 0.5, 40.5, 1000, -0.1, 0.1);


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

	double nentries = 10;
	//double nentries = T->GetEntries();
	for (unsigned int ientry = 0; ientry < nentries; ++ientry) {
		//T->GetEntry(atoi(argv[1]));
		if(ientry%1000==0) std::cout<<"Processing: "<<100.*ientry/nentries <<"%"<<"\n";

		T->GetEntry(ientry);

		if (nhits < 0) {
			//LogDEBUG;
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
		//SMART(genfit::AbsTrackRep) rep = NEW(genfit::RKTrackRep)(pid);
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
		fitter->processTrack(track, false);

		//!
		genfit::MeasuredStateOnPlane* state_at_beam_line = track->extrapolateToLine(
				TVector3(0, 0, 0), TVector3(0, 0, 1));
		//state_at_beam_line->Print();

//		genfit::MeasuredStateOnPlane* state_at_layer_6 = track->extrapolateToCylinder(80.,
//						TVector3(0, 0, 0), TVector3(0, 0, 1));
//		//state_at_layer_6->Print();
//		delete state_at_layer_6;


		TVector3 GenFit_mom = state_at_beam_line->getMom();

		TVector3 AlanDion_mom(AlanDion_px,AlanDion_py,AlanDion_pz);

		hpT_residual_vs_pT->Fill(True_mom.Pt(),(GenFit_mom.Pt() - True_mom.Pt())/True_mom.Pt());

		hDCAr_vs_pT->Fill(True_mom.Mag(),state_at_beam_line->getState()[3]);

		delete state_at_beam_line;
		delete track;
		measurements.clear();
	}

	gStyle->SetOptFit();
	gStyle->SetOptStat(000000000);

	TF1 *tf_pT_resolution = new TF1("tf_pT_resolution","sqrt([0]*[0] + x*x*[1]*[1])", 0, 40);
	tf_pT_resolution->SetParameters(0,0);
	TCanvas *c3 = new TCanvas("c3","c3");
	c3->Divide(2,1);
	c3->cd(1);
	hpT_residual_vs_pT->FitSlicesY();
	TH1D *hpT_resolution_vs_pT = (TH1D*)gDirectory->Get("hpT_residual_vs_pT_2");
	hpT_resolution_vs_pT->SetTitle("PHGenFit: #sigma_{p_{T}}/p_{T}; p_{T}[GeV/c]; #sigma_{p_{T}}/p_{T}");
	hpT_resolution_vs_pT->SetMarkerStyle(20);
	hpT_resolution_vs_pT->Draw("e");
	hpT_resolution_vs_pT->Fit(tf_pT_resolution);
	c3->cd(2);
	hDCAr_vs_pT->FitSlicesY();
	TH1D *hDCAr_resolution_vs_pT = (TH1D*) gDirectory->Get("hDCAr_vs_pT_2");
	hDCAr_resolution_vs_pT->SetTitle(
			"PHGenFit: #sigma_{DCAr} [cm]; p [GeV/c]; #sigma_{DCAr}");
	hDCAr_resolution_vs_pT->SetMarkerStyle(20);
	hDCAr_resolution_vs_pT->Draw("e");
	//hDCAr_resolution_vs_pT->Fit(tf_pT_resolution);
	c3->Print("pT_DCA_resolution.root");

	//! Event display
	//fitter->displayEvent();

	fPHG4Hits->Close();

	delete fitter;

	//pause();

	std::cout<<"SUCCESS! \n";

	return 0;
}
