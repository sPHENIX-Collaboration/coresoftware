//
// Do a recalibration from the saved histograms
//
#include "get_runstr.h"

// Two gaussians
// [0] = ampl, peak 1
// [1] = mean
// [2] = sigma
// [2] = ampl, peak 2
Double_t gaus2(Double_t *x, Double_t *par)
{
   Double_t xx =x[0];
   Double_t f = par[0]*TMath::Gaus(xx,par[1],par[2]) + par[3]*TMath::Gaus(xx,2.0*par[1],sqrt(2)*par[2]); 
   return f;
}

Double_t landau2(Double_t *x, Double_t *par)
{
   Double_t xx =x[0];
   Double_t f = par[0]*TMath::Landau(xx,par[1],par[2]) + par[3]*TMath::Landau(xx,2.0*par[1],par[4]); 
   return f;
}

Double_t landaugaus(Double_t *x, Double_t *par)
{
   Double_t xx =x[0];
   Double_t f = par[0]*TMath::Landau(xx,par[1],par[2]) + par[3]*TMath::Gaus(xx,2.0*par[1],par[4]); 
   return f;
}


void recal_mbd_mip(const char *tfname = "calmbdq_pass1.root", const int pass = 1, const int nevt = 0, const int is_sim = 0)
{
  cout << "tfname " << tfname << endl;

  const int NUM_PMT = 128;
  //const int NUM_PMT = 12;
  const int NUM_ARMS = 2;

  // Read in TFile with h_q
  TFile *oldfile = new TFile(tfname,"READ");

  TH1 *h_q[NUM_PMT];
  TH1 *h_tq[NUM_PMT];

  TString name;
  TString title;
  for (int ipmt=0; ipmt<NUM_PMT; ipmt++)
  {
    name = "h_q"; name += ipmt;
    title = "q"; title += ipmt;
    //h_q[ipmt] = new TH1F(name,title,15100/4,-100,15000);
    h_q[ipmt] = (TH1*)oldfile->Get(name);

    name = "h_tq"; name += ipmt;
    title = "tq"; title += ipmt;
    //h_tq[ipmt] = new TH1F(name,title,7000,-150,31*17.76);
    h_tq[ipmt] = (TH1*)oldfile->Get(name);
  }
  //TH2 *h2_tq = new TH2F("h2_tq","ch vs tq",900,-150,150,NUM_PMT,-0.5,NUM_PMT-0.5);
  TH2 *h2_tq = (TH2*)oldfile->Get("h2_tq");

  // Create new TFile
  TString dir = "results/";
  dir += get_runnumber(tfname);
  dir += "/";
  name = "mkdir -p "; name += dir;
  gSystem->Exec( name );
  name = dir; 
  name += "recalmbd_pass"; name += pass; name += ".root";
  cout << name << endl;

  TFile *savefile = new TFile(name,"RECREATE");

  // Load in calib constants
 
  TCanvas *ac[100];
  int cvindex = 0;

  //tmax
  /*
  ac[cvindex] = new TCanvas("cal_tq","ch vs tq",425*1.5,550*1.5);
  if ( h2_tq != nullptr )
  {
    h2_tq->Draw("colz");
    if ( pass==0 )
    {
      name = dir + "h2_tq.png";
    }
    else
    {
      name = dir + "h2_tqcorr.png";
      h2_tq->GetXaxis()->SetRangeUser(-20,20);
      gPad->Modified();
      gPad->Update();
    }
  }
  cout << name << endl;
  ac[cvindex]->Print( name );
  ++cvindex;
  */

  //q
  ac[cvindex] = new TCanvas("cal_q","q",425*1.5,550*1.5);
  if ( pass>0 )
  {
    ac[cvindex]->Divide(1,2);
    ac[cvindex]->cd(1);
  }

  ofstream cal_mip_file;
  if ( pass==1 ) 
  {
    name = dir + "mbd_qfit.calib";
    cal_mip_file.open( name );
  }

  Double_t qmin = 25.;
  Double_t qmax = 1600;
  if ( is_sim )
  {
    qmin = 0.25;
    qmax = 6.;
  }

  //TF1 *mipfit = new TF1("mipfit","gaus+expo(3)",100,600);
  //TF1 *mipfit = new TF1("mipfit","gaus+pol4(3)",100,600);
  //TF1 *mipfit = new TF1("mipfit","gaus+pol4(3)+expo(8)",qmin,600);
  //TF1 *mipfit = new TF1("mipfit","gaus+pol2(3)+expo(6)",qmin,600);
  //TF1 *mipfit = new TF1("mipfit","gaus+pol3(3)+expo(7)",qmin,600);
  //TF1 *mipfit = new TF1("mipfit","gaussian",qmin,qmax);
  TF1 *mipfit[128];

  /*
     TF1 *f_expo = new TF1("f_expo","expo",qmin,50);
     f_expo->SetParameters(1,-1);
     f_expo->SetLineColor( 2 );

     TF1 *f_bkg = new TF1("f_bkg","pol4+expo(5)",qmin,600);
     f_expo->SetLineColor( 5 );
     */

  TH1 *h_bkg[NUM_PMT];  // background histogram
  TH1 *h_mip[NUM_PMT];  // mip signal histogram
  TH1 *h_bkgmip[NUM_PMT];  // bkg + fit histogram

  //TSpectrum *s = new TSpectrum(1);
  for (int ipmt=0; ipmt<NUM_PMT; ipmt++)
  {

    if (pass>0)
    {
      h_bkg[ipmt] = (TH1*)h_q[ipmt]->Clone();
      name = h_q[ipmt]->GetName(); name.ReplaceAll("q","bkg");
      h_bkg[ipmt]->SetName( name );
      h_bkg[ipmt]->SetTitle( name );
      h_bkg[ipmt]->SetLineColor(2);
      h_q[ipmt]->GetXaxis()->SetRangeUser( qmin, qmax );
      h_q[ipmt]->Sumw2();

      double sigma = 8;
      TSpectrum s{};
      h_bkg[ipmt] = s.Background( h_q[ipmt] );

      h_mip[ipmt] = (TH1*)h_q[ipmt]->Clone();
      name = h_q[ipmt]->GetName(); name.ReplaceAll("q","mip");
      h_mip[ipmt]->SetName( name );
      h_mip[ipmt]->SetTitle( name );
      h_mip[ipmt]->Add( h_bkg[ipmt], -1.0 );

      /*
         Int_t nfound = s.Search(h_q[ipmt],sigma,"",0.1);
         Double_t *xpeaks = s.GetPositionX();

      //h_bkg[ipmt] = s.Background( h_q[ipmt] );

      double best_peak = xpeaks[0];
      if ( best_peak < 50. )
      {
      best_peak = xpeaks[1];
      }

      cout << "peaks\t" << ipmt << "\t" << nfound << "\t" << best_peak << endl;
      */


      // Fit the mip peak after background subtraction
      name = "mipfit"; name += ipmt;
      /*
      // Laudau fit
      mipfit[ipmt] = new TF1(name,"landau",qmin,qmax);
      mipfit[ipmt]->SetLineColor(4);
      mipfit[ipmt]->SetParameter( 0, 5.0*h_mip[ipmt]->GetMaximum() );
      Double_t seedmean = h_mip[ipmt]->GetBinCenter( h_mip[ipmt]->GetMaximumBin() );
      cout << "SEEDMEAN " << seedmean << endl;
      mipfit[ipmt]->SetParameter( 1, seedmean );
      mipfit[ipmt]->SetParameter( 2, 20. );
      */
      // Two gaussian fit
      mipfit[ipmt] = new TF1(name,gaus2,qmin,qmax,4);
      //mipfit[ipmt] = new TF1(name,landau2,qmin,qmax,5);
      //mipfit[ipmt] = new TF1(name,landaugaus,qmin,qmax,5);
      mipfit[ipmt]->SetLineColor(4);
      mipfit[ipmt]->SetParameter( 0, 5.0*h_mip[ipmt]->GetMaximum() );
      Double_t seedmean = h_mip[ipmt]->GetBinCenter( h_mip[ipmt]->GetMaximumBin() );
      cout << "SEEDMEAN " << seedmean << endl;
      Double_t seedsigma = 20.;
      if ( is_sim ) seedsigma = 0.2;
      mipfit[ipmt]->SetParameter( 1, seedmean );
      mipfit[ipmt]->SetParameter( 2, seedsigma );
      mipfit[ipmt]->SetParameter( 3, mipfit[ipmt]->GetParameter(0)*0.1 );
      //mipfit[ipmt]->SetParameter( 4, mipfit[ipmt]->GetParameter(1) );

      h_mip[ipmt]->Fit( mipfit[ipmt], "RM" );

      double integ = mipfit[ipmt]->GetParameter(0);
      double best_peak = mipfit[ipmt]->GetParameter(1);
      double width = mipfit[ipmt]->GetParameter(2);
      double integerr = mipfit[ipmt]->GetParError(0);
      double best_peakerr = mipfit[ipmt]->GetParError(1);
      double widtherr = mipfit[ipmt]->GetParError(2);
      double chi2 = mipfit[ipmt]->GetChisquare();
      double ndf = mipfit[ipmt]->GetNDF();

      cal_mip_file << ipmt << "\t" << integ << "\t" << best_peak << "\t" << width << "\t"
        << integerr << "\t" << best_peakerr << "\t" << widtherr << "\t"
        << chi2/ndf << endl;
      cout << ipmt << "\t" << integ << "\t" << best_peak << "\t" << width << "\t"
        << integerr << "\t" << best_peakerr << "\t" << widtherr << "\t"
        << chi2/ndf << endl;

      // Get full fit
      h_bkgmip[ipmt] = (TH1*)h_bkg[ipmt]->Clone();
      name = h_q[ipmt]->GetName(); name.ReplaceAll("q","bkgmip");
      h_bkgmip[ipmt]->SetName( name );
      h_bkgmip[ipmt]->SetTitle( name );
      h_bkgmip[ipmt]->Add( mipfit[ipmt] );
      h_bkgmip[ipmt]->SetLineColor( 8 );

      /*
         gPad->Modified();
         gPad->Update();
      //mipfit->SetRange( xpeaks[0]-2.5*sigma, 600 );
      mipfit->SetRange( qmin, 600 );
      mipfit->SetParameter( 0, 1000 );
      mipfit->SetParameter( 1, best_peak );
      mipfit->SetParameter( 2, sigma );

      mipfit->SetParameter( 8, f_expo->GetParameter(0) );
      mipfit->SetParameter( 9, f_expo->GetParameter(1) );
      h_q[ipmt]->Fit( mipfit, "R" );
      */


      /*
         int npar = f_bkg->GetNpar();
         for (int ipar=0; ipar<npar; ipar++)
         {
         f_bkg->SetParameter( ipar, mipfit->GetParameter(ipar+3) );
         }
         */

      // Now draw the full dist, plus fit
      ac[cvindex]->cd(1);
      gPad->SetLogy(1);
      h_q[ipmt]->GetXaxis()->SetRangeUser( qmin, qmax );
      //h_q[ipmt]->SetMinimum(10.);
      h_q[ipmt]->Draw();
      h_bkg[ipmt]->Draw("same");
      h_bkgmip[ipmt]->Draw("same");

      gPad->Modified();
      gPad->Update();

      ac[cvindex]->cd(2);
      h_mip[ipmt]->Draw();
      gPad->Modified();
      gPad->Update();

      /*
      string junk;
      cout << "? ";
      cin >> junk;
      */

      name = dir + "h_qfit"; name += ipmt; name += ".png";
      cout << name << endl;
      ac[cvindex]->Print( name );
    }

  }
  ++cvindex;

  if ( pass==1 )
  {
    cal_mip_file.close();
  }

  savefile->Write();
  savefile->Close();
}

