/*
   Authors: Haiwang Yu
*/
#include "Field2D.h"

#include <TFile.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH2D.h>

#include <iostream>
#include <cmath>

#define LogDEBUG std::cout<<"DEBUG: "<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"\n"

//int round (double const x){
//	if(x>=0){
//		return (int) (x+0.5);
//	} else {
//		return (int) (x-0.5);
//	}
//
//}

namespace genfit {

bool Field2D::re_scale(double r)
{
	field_map_r_->Scale(r);
	field_map_z_->Scale(r);

	return true;
}

bool Field2D::initialize(std::string inname) {


	field_map_r_ = new TH2D("field_map_r_","B_{r} [kGauss]; z [cm]; r [cm]", 401, -401, 401, 151, -1, 301);
	field_map_z_ = new TH2D("field_map_z_","B_{z} [kGauss]; z [cm]; r [cm]", 401, -401, 401, 151, -1, 301);
	TFile *fin = TFile::Open(inname.data(),"READ");
	if(!fin)
	{
		//TODO Output some msg
		LogDEBUG;
		return false;
	}
	TNtuple *T = (TNtuple *) fin->Get("fieldmap");
	if(!T)
	{
		//TODO Output some msg
		LogDEBUG;
		return false;
	}

	Float_t r;
	Float_t z;
	Float_t br;
	Float_t bz;
	T->SetBranchAddress("r",&r);
	T->SetBranchAddress("z",&z);
	T->SetBranchAddress("br",&br);
	T->SetBranchAddress("bz",&bz);

	for(long ientry=0;ientry<T->GetEntries();ientry++){
		//if (T->GetEntry(ientry) < 0) continue;
		T->GetEntry(ientry);
		//std::cout<<r <<" ,"<<z <<" ,"<<br<<" ,"<<bz<<"\n";
		field_map_r_->SetBinContent(
				field_map_r_->GetXaxis()->FindBin(z),
				field_map_r_->GetYaxis()->FindBin(r),
				br*10);
		field_map_z_->SetBinContent(
				field_map_z_->GetXaxis()->FindBin(z),
				field_map_z_->GetYaxis()->FindBin(r),
				bz*10);
	}
	fin->Close();

	//field_map_z_->Print("all");

	return true;
}

void Field2D::plot(std::string option){

//	TH2D *hbr = new TH2D("hbr","|B_{r}| [kGauss]; z [cm]; r [cm]",401, -401, 401, 151, -1, 301);
//	TH2D *hbz = new TH2D("hbz","|B_{z}| [kGauss]; z [cm]; r [cm]",401, -401, 401, 151, -1, 301);
	TH2D *hbr = new TH2D("hbr","|B_{r}| [kGauss]; z [cm]; r [cm]",300, -350, 350, 200, 0, 290);
	TH2D *hbz = new TH2D("hbz","|B_{z}| [kGauss]; z [cm]; r [cm]",300, -350, 350, 200, 0, 290);
	for (double r = 0; r < 300; r+=1)
		for (double z = 400; z > -400; z-=1) {
			double bx;
			double by;
			double bz;
			get(r, 0, z, bx, by, bz);

//			double br = sqrt(bx * bx + by * by);

			//std::cout <<"DEBUG: "<<__LINE__<<":"<< r << " ," << z << " ," << br << " ," << bz << "\n";

			hbr->SetBinContent(
					hbr->GetXaxis()->FindBin(z),
					hbr->GetYaxis()->FindBin(r),
					bz*10);
			hbz->SetBinContent(
					hbz->GetXaxis()->FindBin(z),
					hbz->GetYaxis()->FindBin(r),
					bz*10);
		}

	TCanvas *c0 = new TCanvas("c0","c0");
	c0->Divide(1,2);
	c0->cd(1);
	hbr->SetStats(0);
	hbr->Draw("colz");
	c0->cd(2);
	hbz->SetStats(0);
	hbz->Draw("colz");
	c0->Update();

	c0->SaveAs("Field2D_plot.root");
	c0->SaveAs("Field2D_plot.pdf");

	delete c0;
	delete hbr;
	delete hbz;

	return;
}

TVector3 Field2D::get(const TVector3& v) const {
  double x = v.x();
  double y = v.y();
  double z = v.z();
  double Bx;
  double By;
  double Bz;
  get(x,y,z,Bx,By,Bz);
  return TVector3(Bx,By,Bz);
}

void Field2D::get(const double&x, const double&y, const double&z, double& Bx, double& By, double& Bz) const {
	double r = sqrt(x*x + y*y);

	if(fabs(r) > 300 || fabs(z) > 400)
	{
		Bx = 0;
		By = 0;
		Bz = 0;
		return;
	}

	int bin_z = field_map_r_->GetXaxis()->FindBin(z);
	int bin_r = field_map_r_->GetYaxis()->FindBin(r);

	Bz = field_map_z_->GetBinContent(bin_z,bin_r);
	double Br = field_map_r_->GetBinContent(bin_z,bin_r);

	Bx = x/r*Br;
	By = y/r*Br;

	if(r == 0)
	{
		Bx = 0;
		By = 0;
	}



	//std::cout<<"DEBUG: "<<__LINE__<<": "<<z<<","<<r<<","<<bin_z<<","<<bin_r<<": "<<Br<<","<<Bz<<"\n";
}

} /* End of namespace genfit */
