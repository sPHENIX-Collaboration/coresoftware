#include "TFile.h"
#include "TTree.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"

#include <stdio.h>
#include <iostream>

#include <cdbobjects/CDBTTree.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libcdbobjects.so)

std::vector<float> vec_peak_pos_val;
std::vector<float> vec_peak_pos_err;

int fit_dz(const char* filename = "data.root",
    const char* treename = "tree",
    const char* branchname = "dz",
    TString cut = "1",
    TString outfile = "fit_dz.pdf",
    TString title = "Run 51103, fit dz",
    int runnumber = 51103)
{
    TChain* chain = new TChain(treename);
    chain->Add(filename);
    int nevent = chain->GetEntries();

    if (nevent > 0) {
        std::cout << "Successfully added file to the TChain." << std::endl;
    } else {
        std::cout << "Failed to add file to the TChain." << std::endl;
        return -1;
    }

    if (chain->GetEntries(cut) < 50)
    {
      std::cout<<"Cut: "<<cut<<" , nevent used in fit = "<<chain->GetEntries(cut)<<" . Too Low!!! Fit performance maybe not be good, skipping!!!"<<std::endl;
      return -1;
    }

    double xmin = -20;
    double xmax = 20;
    int nbins = 50;
    RooRealVar x(branchname, branchname, nbins, xmin, xmax);
    double binwidth = (xmax-xmin)/nbins;

    TFile* newfile = new TFile("tmp.root","recreate");
    TTree* newtree = (TTree*) chain->CopyTree(cut);

    // Create a RooDataSet from the TTree
    RooDataSet data("data", "dataset of dz", newtree, x);

    // Define the parameters for the Gaussian
    RooRealVar mean("mean", "mean of gaussian", data.mean(x), xmin, xmax);
    RooRealVar sigma("sigma", "width of gaussian", 0.1*data.sigma(x), 0.001, (xmax-xmin)/2);

    // Create the Gaussian PDF
    RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);
    RooRealVar ngauss("ngauss","number of gaussian",0.2*data.sumEntries(),0,data.sumEntries());

    // Define the parameters for the linear background
    RooRealVar a0("a0", "constant", 0, -10, 10);
    RooRealVar a1("a1", "slope", 0, -1, 1);
    RooPolynomial poly("poly", "linear background", x, RooArgList(a1, a0));
    RooRealVar npoly("npoly","number of poly",0.3*data.sumEntries(),0,data.sumEntries());

    // Combine the Gaussian and the polynomial into a single PDF
    RooAddPdf model("model", "gauss + poly", RooArgList(gauss, poly), RooArgList(ngauss, npoly));

    // Fit the Gaussian to the data
    RooFitResult* fit_result = model.fitTo(data, RooFit::Save());

    // Print the fit result
    fit_result->Print();

    // Plot the data and the fit result
    RooPlot* frame = x.frame();
    data.plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, RooFit::Components(gauss), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
    model.plotOn(frame, RooFit::Components(poly), RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));

    frame->SetTitle(title);
    frame->SetXTitle("#DeltaZ [cm]");
    frame->SetYTitle(Form("Events / (%1.0g cm)",binwidth));
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);

    double meanValue = mean.getVal();
    double meanError = mean.getError();
    double sigmaValue = sigma.getVal();
    double sigmaError = sigma.getError();

    vec_peak_pos_val.push_back(meanValue);
    vec_peak_pos_err.push_back(meanError);

/*
    TPaveText *pt = new TPaveText(0.20, 0.73, 0.35, 0.88, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->AddText(Form("#mu = %.2f #pm %.2f", meanValue, meanError));
    pt->AddText(Form("#sigma = %.2f #pm %.2f", sigmaValue, sigmaError));

    // Draw the plot on a canvas
    TCanvas *canvas = new TCanvas("canvas", "fit result", 800, 600);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.05);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.15);
    frame->Draw();
    pt->Draw();

    std::string dirName = "figure/" + std::to_string(runnumber);
    std::string command = "ls " + dirName + " 2>/dev/null";
    int result = system(command.c_str());

    if (result != 0) {
        std::string mkdirCommand = "mkdir -p " + std::string(dirName);
        int mkdirResult = system(mkdirCommand.c_str());

        if (mkdirResult == 0) {
            std::cout << "Directory '" << dirName << "' created successfully." << std::endl;
        } else {
            std::cerr << "Failed to create directory '" << dirName << "'." << std::endl;
        }
    } else {
        std::cout << "Directory '" << dirName << "' already exists." << std::endl;
    }

    canvas->SaveAs(outfile);
    delete canvas;
*/

    delete chain;

    system("rm tmp.root");
    return 1;
}

void CDBTTreeWrite(const std::string &fname = "test.root", float dv=0.007)
{
  CDBTTree *cdbttree = new CDBTTree(fname);
  cdbttree->SetSingleFloatValue("tpc_drift_velocity",dv); // unit: cm/ns
  cdbttree->CommitSingle();
  cdbttree->Print();
  cdbttree->WriteCDBTTree();
  delete cdbttree;
  //gSystem->Exit(0);
}

void fit(int runnumber)
{
    gStyle->SetOptStat(0);
    std::vector<float> vec_trkr_z_val;
    std::vector<float> vec_trkr_z_err;

    // fit dz in different trkr_z region
    vec_peak_pos_val.clear();
    vec_peak_pos_err.clear();
    vec_trkr_z_val.clear();
    vec_trkr_z_err.clear();
    float z_min = -140;
    float z_max = 140;
    float z_range = z_max - z_min;
    int nstep = 14;
    float z_stepsize = z_range / nstep;
    for (int i=0; i<nstep; i++)
    {
      float z_cut_min = z_min + i * z_stepsize;
      float z_cut_max = z_min + (i+1) * z_stepsize;
      int status = fit_dz(Form("Reconstructed/%d/TrackCalo_*_ana.root",runnumber), "tree", "dz", Form("fabs(dphi)<0.1 && track_z>%f && track_z<%f && fabs(dz)<20",z_cut_min,z_cut_max), Form("figure/%d/dz_fit_cuttrkrz_%d.pdf",runnumber,i), Form("Run %d, Track Z#in[%.0f,%.0f] cm",runnumber,z_cut_min,z_cut_max),runnumber);
      if (status==1)
      {
        vec_trkr_z_val.push_back((z_cut_min+z_cut_max)/2);
        vec_trkr_z_err.push_back((z_cut_max-z_cut_min)/2);
      }
    }

    int npoint = vec_peak_pos_val.size();

    if (npoint==0)
    {
      std::cout<<"Do not have anything to fit!!! DST is empty or Matching for this run is very bad"<<std::endl;
      return;
    }

    int npoint_left=0;
    int npoint_right=0;
    for (int i=0; i<(vec_peak_pos_val.size()); i++)
    {
      if (vec_trkr_z_val.at(i)>0)
      {
        npoint_right++;
      }
      else if (vec_trkr_z_val.at(i)<0)
      {
        npoint_left++;
      }
    }
    if (npoint_left<3 || npoint_right<3)
    {
      std::cout<<"Do not have enough good point to fit!!! Matching for this run is very bad"<<std::endl;
      return;
    }

    TCanvas *c1 = new TCanvas("c1", "Fitting Example", 800, 600);
    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.05);
    TGraphErrors *graph = new TGraphErrors(vec_peak_pos_val.size(), vec_trkr_z_val.data(), vec_peak_pos_val.data(), vec_trkr_z_err.data(), vec_peak_pos_err.data());
    graph->SetTitle(Form("Run %d, TPC-EMCal matching;Track Z [cm];#DeltaZ=Z_{track}-Z_{calo} Peak [cm]",runnumber));
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    graph->GetYaxis()->SetRangeUser(-50,20);

    double intercept_min = -20;
    double intercept_max= 20;
    double slope_min = -1;
    double slope_max = 1;

    TF1 *fitFunc1 = new TF1("fitFunc1", "[0] + [1]*x", z_min, -10);
    fitFunc1->SetParameters(vec_peak_pos_val.at(npoint/2 - 2), (vec_peak_pos_val.at(npoint/2 - 2) - vec_peak_pos_val.at(0)) / (vec_trkr_z_val.at(npoint/2 - 2) - vec_trkr_z_val.at(0)));
    fitFunc1->SetParLimits(0, intercept_min, intercept_max);
    fitFunc1->SetParLimits(1, slope_min, slope_max);
    graph->Fit(fitFunc1, "R");

    TF1 *fitFunc2 = new TF1("fitFunc2", "[0] + [1]*x", 10, z_max);
    fitFunc2->SetParameters(vec_peak_pos_val.at(npoint/2 + 1), (vec_peak_pos_val.at(npoint - 1) - vec_peak_pos_val.at(npoint/2 + 1)) / (vec_trkr_z_val.at(npoint - 1) - vec_trkr_z_val.at(npoint/2 + 1)));
    fitFunc2->SetParLimits(0, intercept_min, intercept_max);
    fitFunc2->SetParLimits(1, slope_min, slope_max);
    graph->Fit(fitFunc2, "R");

    double chi2_1 = fitFunc1->GetChisquare();
    int ndf_1 = fitFunc1->GetNDF();

    double chi2_2 = fitFunc2->GetChisquare();
    int ndf_2 = fitFunc2->GetNDF();

    double p0_1 = fitFunc1->GetParameter(0);
    double p1_1 = fitFunc1->GetParameter(1);
    double p0_1_err = fitFunc1->GetParError(0);
    double p1_1_err = fitFunc1->GetParError(1);
    TString txt_fun1 = Form("y = %.3fx + (%.3f)",p1_1,p0_1);
    TString txt_quality1 = Form("#chi^2/NDF = %.2f/%d = %.2f",chi2_1,ndf_1,chi2_1/ndf_1);
    std::cout << "Fit left part: " << txt_fun1 << " , " << txt_quality1 << std::endl;

    double p0_2 = fitFunc2->GetParameter(0);
    double p1_2 = fitFunc2->GetParameter(1);
    double p0_2_err = fitFunc2->GetParError(0);
    double p1_2_err = fitFunc2->GetParError(1);
    TString txt_fun2 = Form("y = %.3fx + (%.3f)",p1_2,p0_2);
    TString txt_quality2 = Form("#chi^2/NDF = %.2f/%d = %.2f",chi2_2,ndf_2,chi2_2/ndf_2);
    std::cout << "Fit right part: " << txt_fun2 << " , " << txt_quality2 << std::endl;

    graph->Draw("AP");

    fitFunc1->SetLineColor(kRed);
    fitFunc1->Draw("SAME");
    fitFunc2->SetLineColor(kBlue);
    fitFunc2->Draw("SAME");

    double driftvelo_pre = 0.00710;
    double dz_separation_0_new = -(p1_1+p1_2)/2. * 2 * 105;
    double dz_separation_0_err = -(p1_1_err+p1_2_err)/2. * 2 * 105;
    double driftvelo_new = driftvelo_pre * (1+ dz_separation_0_new / 2. / 105.);
    double driftvelo_err = driftvelo_pre * (fabs(dz_separation_0_err) / 2. / 105.);
    TString printtxt1 = Form("DV used in reconstruction: %.5f cm/ns",driftvelo_pre);
    TString printtxt2 = Form("Calibrated DV from data: (%.5f#pm%.5f) cm/ns",driftvelo_new,driftvelo_err);
    cout<<printtxt1<<endl;
    cout<<printtxt2<<endl;

    TPaveText *pt1 = new TPaveText(0.15, 0.275, 0.45, 0.40, "NDC");
    pt1->SetFillColor(0);
    pt1->SetTextAlign(12);
    pt1->AddText(txt_fun1);
    pt1->AddText(txt_quality1);
    pt1->Draw("same");

    TPaveText *pt2 = new TPaveText(0.55, 0.275, 0.85, 0.40, "NDC");
    pt2->SetFillColor(0);
    pt2->SetTextAlign(12);
    pt2->AddText(txt_fun2);
    pt2->AddText(txt_quality2);
    pt2->Draw("same");

    TPaveText *pt3 = new TPaveText(0.20, 0.15, 0.90, 0.275, "NDC");
    pt3->SetFillColor(0);
    pt3->SetTextAlign(12);
    pt3->AddText(printtxt1);
    pt3->AddText(printtxt2);
    pt3->Draw("same");

    TPaveText *pt4 = new TPaveText(0.15, 0.40, 0.35, 0.50, "NDC");
    pt4->SetFillColor(0);
    pt4->SetTextAlign(12);

    // checks if result is good to insert into CDB
    double fit_quality_thr = 10;
    bool fit_status = true;
    if (chi2_1/ndf_1>fit_quality_thr || chi2_2/ndf_2>fit_quality_thr)
    {
      cout<<"ERROR: fail fit quality cut"<<endl;
      fit_status = false;
    }

    if (fabs(p0_1-intercept_min)<1e-3 || fabs(p0_1-intercept_max)<1e-3
     || fabs(p1_1-slope_min)<1e-3 || fabs(p1_1-slope_max)<1e-3
     || fabs(p0_2-intercept_min)<1e-3 || fabs(p0_2-intercept_max)<1e-3
     || fabs(p1_2-slope_min)<1e-3 || fabs(p1_2-slope_max)<1e-3)
    {
      cout<<"ERROR: parameters in the limit"<<endl;
      fit_status = false;
    }

    if ( fabs(p1_1 - p1_2) > 0.1 )
    {
      cout<<"ERROR: left fit and right fit inconsistent"<<endl;
      fit_status = false;
    }

    double dv_ref = driftvelo_pre;
    if (runnumber > 50800 && runnumber < 51200)
    {
      dv_ref = 0.00625; // for TPC gas mixture ratio changes during 08/09 - 08/12
    }
    if (driftvelo_new/dv_ref > 1.1 || driftvelo_new/dv_ref < 0.9)
    {
      fit_status = false;
    }

    if (fit_status)
    {
      pt4->SetTextColor(3);
      pt4->AddText("GOOD");
      pt4->Draw("same");
      c1->Update();
      c1->SaveAs(Form("./figure/dz_fit_trkrz_Run%d.pdf",runnumber));
      CDBTTreeWrite(Form("cdbttree/tpc_drift_velocity_%d.root",runnumber), (float) driftvelo_new);
    }
    else
    {
      pt4->SetTextColor(2);
      pt4->AddText("BAD");
      pt4->Draw("same");
      c1->Update();
      c1->SaveAs(Form("./figure_failed/dz_fit_trkrz_Run%d.pdf",runnumber));
      CDBTTreeWrite(Form("cdbttree_failed/tpc_drift_velocity_%d.root",runnumber), (float) driftvelo_new);
    }
}
