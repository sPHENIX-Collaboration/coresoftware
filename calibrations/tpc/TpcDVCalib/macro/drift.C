/*
 * This macro shows a minimum working example of running the tracking
 * hit unpackers with some basic seeding algorithms to try to put together
 * tracks. There are some analysis modules run at the end which package
 * hits, clusters, and clusters on tracks into trees for analysis.
 */

#include <fun4all/Fun4AllUtils.h>
// #include <G4_ActsGeom.C>
// #include <G4_Global.C>
// #include <G4_Magnet.C>
// #include <G4_Mbd.C>
#include <GlobalVariables.C>
// #include <QA.C>
// #include <Trkr_Clustering.C>
// #include <Trkr_Reco.C>
// #include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
// #include <fun4all/Fun4AllDstInputManager.h>
// #include <fun4all/Fun4AllDstOutputManager.h>
// #include <fun4all/Fun4AllInputManager.h>
// #include <fun4all/Fun4AllOutputManager.h>
// #include <fun4all/Fun4AllRunNodeInputManager.h>
// #include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

// #include <trackingqa/InttClusterQA.h>
// #include <trackingqa/MicromegasClusterQA.h>
// #include <trackingqa/MvtxClusterQA.h>
// #include <trackingqa/TpcClusterQA.h>
// 
// #include <trackingdiagnostics/TrackResiduals.h>
// #include <trackingdiagnostics/TrkrNtuplizer.h>

#include <stdio.h>

// #include <collect/Collect.h>

//R__LOAD_LIBRARY(libfun4all.so)
//R__LOAD_LIBRARY(libffamodules.so)
//R__LOAD_LIBRARY(libmvtx.so)
//R__LOAD_LIBRARY(libintt.so)
R__LOAD_LIBRARY(libtpc.so)
//R__LOAD_LIBRARY(libmicromegas.so)
//R__LOAD_LIBRARY(libTrackingDiagnostics.so)
//R__LOAD_LIBRARY(libtrackingqa.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
//R__LOAD_LIBRARY(libcollect.so)

double driftGET( const int arun = 51279 );
bool CheckHV( const int arun = 51279 );

using std::cout;
using std::endl;

//------------------------------------------------------------------------
//	get drift velocity from cdb...	wjl aug 2024
//
void drift(){

//	~70 um/ns -> 0.007 in sphenix units

	TGraph *gdriftrun	= new TGraph();
			gdriftrun	->SetMarkerStyle(20);
			gdriftrun	->SetMarkerSize(0.9);
			gdriftrun	->SetMarkerColor(kBlue);
			gdriftrun	->SetLineColor(kBlue);

	double drifts[2000]	= {0};
	const int NPHASE = 9;
	int iphase = 0;
	const char* phasename[NPHASE]	= {"6x6","28x28","56x56","111x111 8/9","111x111 8/11","","BP#downarrow","",""};
	int    lastruns[NPHASE]	= {49722,50025,50458,50936,51109,51191,51301,51495,51617};
	double phaseN[NPHASE]	= {0};
	double phaseA[NPHASE]	= {0};
	double phaseNhv[NPHASE]	= {0};
	double phaseAhv[NPHASE]	= {0};
	double lowest[NPHASE]	= {0};
	for (int iph=0;iph<NPHASE;iph++){ lowest[iph] = 99.; }

	//---- loop over run numbers...
	//
	for (int i=0;i<2000;i++){
		//
		int irun	= 49700+i;
		double val	= driftGET(irun);
		//
		if (irun<=49722){ iphase=0; } else
		if (irun<=50025){ iphase=1; } else
		if (irun<=50458){ iphase=2; } else
		if (irun<=50936){ iphase=3; } else
		if (irun<=51109){ iphase=4; } else
		if (irun<=51191){ iphase=5; } else
		if (irun<=51301){ iphase=6; } else 
		if (irun<=51495){ iphase=7; } else 
                        { iphase=8; } 
        //
        bool CorrectHV	= CheckHV(irun);
        //
		drifts[i]	= val;
		if (val>0.){
			gdriftrun	->SetPoint(gdriftrun->GetN(),irun,val);
			phaseA[iphase]	+= val;
			phaseN[iphase]	+= 1.0;
//			if (CorrectHV){
//				phaseAhv[iphase]	+= val;
//				phaseNhv[iphase]	+= 1.0;
//			}
			if (val<lowest[iphase]) lowest[iphase]=val;
		}
	}
	
	int nrunstot	= 0;
	int nrunstothv	= 0;
	for (int iph=0;iph<NPHASE;iph++){
		if (phaseN[iph]>0){
			nrunstot		+= phaseN[iph];
			nrunstothv		+= phaseNhv[iph];
			phaseA[iph]		/= phaseN[iph];
//			phaseAhv[iph]	/= phaseNhv[iph];
//			cout<<"RunNum<="<<lastruns[iph]<<"\t Nruns="<<phaseN[iph]<<" "<<phaseNhv[iph]<<"\t <v>="<<phaseA[iph]<<endl;
		}
	}
	cout<<"Nruns with custom <v> found = "<<nrunstot<<endl;
		
	//---- paint setup
	int ican=-1,itext=-1,iline=-1,ivline=-1,iframe=-1; 
	TCanvas	*ccan[100];
	TH1F	*frame[100];
	TLatex 	*text[1000];
	TLine 	*line[100];
	TLine 	*vline[100];
	for (int i=0;i<1000;i++){
		text[i]	= new TLatex();
		text[i]	->SetNDC();
		text[i]	->SetTextFont(42);
		text[i]	->SetTextSize(0.05);
		text[i]	->SetTextAlign(33);
		line[i]	= new TLine();
		line[i]	->SetLineColor(1);
		line[i]	->SetLineStyle(2);
		line[i]	->SetLineWidth(2);
		vline[i]= new TLine(0,-1,0,1);
		vline[i]->SetLineColor(17);
		vline[i]->SetLineWidth(2);
	}
	gROOT->SetStyle("Modern");
	gStyle->SetOptStat(0);
 	gStyle->SetTitleAlign(13);
	gStyle->SetTitleFontSize(0.08);		
 	gStyle->SetTitleX(0.20);
 	gStyle->SetTitleY(0.98);
	gStyle->SetPadRightMargin(0.01);
	gStyle->SetPadTopMargin(0.01);
	gStyle->SetPadBottomMargin(0.08);
	gStyle->SetPadLeftMargin(0.12);
 	gStyle->SetTitleSize(0.03,"xyzt");
 	gStyle->SetLabelSize(0.03,"xyzt");
 	//
 	//---- end paint setup..

	//---- paint...
	//	
	++ican; ccan[ican]	= new TCanvas(Form("ccan%d",ican),Form("ccan%d",ican),20+ican*30,30+ican*30,0.7*1100,0.7*850);
	ccan[ican]->cd(); ccan[ican]->Divide(1,1,0.0001,0.0001);
		ccan[ican]->cd(1);
		gdriftrun->Draw("AP");
		for (int iph=0;iph<NPHASE;iph++){
			++itext; text[itext]->SetNDC(kFALSE); text[itext]->SetTextAlign(23); text[itext]->SetTextSize(0.03); 
				text[itext]->DrawLatex(lastruns[iph]-20,0.996*lowest[iph],Form("%.5f",phaseA[iph]));
			//++itext; text[itext]->SetNDC(kFALSE); text[itext]->SetTextAlign(23); text[itext]->SetTextSize(0.03); text[itext]->SetTextColor(kGreen+1); 
			//	text[itext]->DrawLatex(lastruns[iph],0.996*lowest[iph]-0.00005,Form("%.5f",phaseAhv[iph]));
			//
			++itext; text[itext]->SetNDC(kFALSE); text[itext]->SetTextAlign(21); text[itext]->SetTextSize(0.03); 
				text[itext]->SetTextAngle(90);
				text[itext]->DrawLatex(lastruns[iph],0.0067,Form("%s",phasename[iph]));
		}
	ccan[ican]->cd(); ccan[ican]->Update();
	ccan[ican]->Print("drift.ps(");
	
	//---- close-out...
	//
	ccan[ican]->Print("drift.ps]");
	int iSuccess = gSystem->Exec("ps2pdf drift.ps drift.pdf");
	cout<<"isuccess="<<iSuccess<<endl;
	if (iSuccess==-1){
		cout<<"drift.pdf produced, deleting drift.ps file..."<<endl;
		TString sexec	 = TString("/bin/rm drift.ps");
		cout<<sexec.Data()<<endl;
		gSystem			->Exec(sexec.Data());
	}
	cout<<"Done..."<<endl;
	
	cout<<"Writing drift.root..."<<endl;
	TFile *fout	= new TFile("drift.root","RECREATE");
		fout->cd();
		gdriftrun->Write("gdriftrun");
	fout->Close();
		
}

//------------------------------------------------------------
//
//	return true if HV=3.3 kV
//
bool CheckHV( const int irun = 51279 ){
	if (irun>=49706&&irun<=49713){ return true; } else
	if (irun>=50006&&irun<=50014){ return true; } else
	if (irun>=50021&&irun<=50025){ return true; } else
	if (irun>=50438&&irun<=50440){ return true; } else
	if (irun>=50444&&irun<=50449){ return true; } else
	if (irun>=50916&&irun<=50921){ return true; } else
	if (irun>=50933&&irun<=50936){ return true; } else
	if (irun>=51099&&irun<=51107){ return true; }
// 	if (irun>=&&irun<=){ return true; } else
// 	if (irun>=&&irun<=){ return true; } else
// 	if (irun>=&&irun<=){ return true; } else
// 	if (irun>=&&irun<=){ return true; } else
	return false;
}

//-----------------------------------------------------------
//
double driftGET( const int runnumber = 51279 ){
	//
	TpcReadoutInit( runnumber );
	double driftOLD	= G4TPC::tpc_drift_velocity_reco;
	//
	auto rc = recoConsts::instance();
	rc->set_IntFlag("RUNNUMBER", runnumber);
	//	
	Enable::CDB = true;
	rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
	rc->set_uint64Flag("TIMESTAMP", 6);
	//
	//G4TPC::tpc_drift_velocity_reco = (8.0 / 1000) * 107.0 / 105.0;
	//
	//---- get drift velocity
	double driftCDB	= 0;
	rc->set_uint64Flag("TIMESTAMP",runnumber);
 	CDBInterface* cdb = CDBInterface::instance();			
	std::string tpc_dv_calib_dir = cdb->getUrl("TPC_DRIFT_VELOCITY");
	if (!tpc_dv_calib_dir.empty()){			// path to driftvelocity calib file was found...
		//cout << "DRIFT: Path of tpc_drift_velocity_calib_file = "<< tpc_dv_calib_dir << endl;
		CDBTTree *cdbttree = new CDBTTree(tpc_dv_calib_dir);
		cdbttree->LoadCalibrations();
		driftCDB	= cdbttree->GetSingleFloatValue("tpc_drift_velocity");
		cout<<"DRIFT: RunNum="<<runnumber<<"\t driftOLD=" << driftOLD << "\t driftCDB=" << driftCDB << endl;
		//if (driftCDB>0.&&driftCDB<0.01){
		//	G4TPC::tpc_drift_velocity_reco	= driftCDB;
		//}
	} else {
		//cout << "DRIFT: Path to tpc_drift_velocity_calib_file NOT FOUND...."<< endl;
	}
	//
	return driftCDB;
}

