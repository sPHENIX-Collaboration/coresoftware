#include <TF1.h>

#include "get_runstr.h"
#include "mbd/MbdCalib.h"
#include "mbd/MbdDefs.h"
#include "read_dstmbd.C"
#include "recal_mbd_mip.C"


#if defined(__CLING__)
R__LOAD_LIBRARY(libmbd.so)
#endif

//pass0: do tt_t0 and tq_t0 offset calibration
//pass1: 
void cal_mbd(const char *tfname = "DST_MBDUNCAL-00020869-0000.root", const int pass = 0, const int nevt = 0)
{
  cout << "tfname " << tfname << endl;
  read_dstmbd( tfname );

  //== Create output directory (should already exist)
  TString dir = "results/";
  dir += get_runnumber(tfname);
  dir += "/";
  TString name = "mkdir -p " + dir;
  gSystem->Exec( name );

  //== Load in calib constants
  //Float_t tq_t0_offsets[NUM_PMT] = {};

  MbdCalib *mcal = new MbdCalib();
  TString calfile = dir + "/mbd_sampmax.calib";
  mcal->Download_SampMax( calfile.Data() );

  TString savefname = dir;
  if ( pass==0 )
  {
    savefname += "calmbdtime_pass"; savefname += pass; savefname += ".root";
  }
  else if ( pass>0 )
  {
    calfile = dir + "/mbd_tq_t0.calib";
    mcal->Download_TQT0( calfile.Data() );

    calfile = dir + "/mbd_tt_t0.calib";
    mcal->Download_TTT0( calfile.Data() );

    savefname += "calmbdq_pass"; savefname += pass; savefname += ".root";
  }
  if ( pass>1 )
  {
    calfile = dir + "/mbd_qfit.calib";
    mcal->Download_Gains( calfile.Data() );

    savefname += "calmbd_pass"; savefname += pass; savefname += ".root";
  }

  TFile *savefile = new TFile(savefname,"RECREATE");

  TH1 *h_q[NUM_PMT];
  TH1 *h_tt[NUM_PMT];
  TH1 *h_tq[NUM_PMT];

  TString title;
  for (int ipmt=0; ipmt<NUM_PMT; ipmt++)
  {
    name = "h_q"; name += ipmt;
    title = "q"; title += ipmt;
    h_q[ipmt] = new TH1F(name,title,3000,-100,15000-100);

    name = "h_tt"; name += ipmt;
    title = "tt"; title += ipmt;
    h_tt[ipmt] = new TH1F(name,title,7000,-150,31*17.76);

    name = "h_tq"; name += ipmt;
    title = "tq"; title += ipmt;
    h_tq[ipmt] = new TH1F(name,title,7000,-150,31*17.76);
  }
  TH2 *h2_tq = new TH2F("h2_tq","ch vs tq",900,-150,150,NUM_PMT,-0.5,NUM_PMT-0.5);
  TH2 *h2_tt = new TH2F("h2_tt","ch vs tt",900,-150,150,NUM_PMT,-0.5,NUM_PMT-0.5);

  // Event loop, each ientry is one triggered event
  int nentries = tree->GetEntries();
  if ( nevt!=0 && nevt<nentries )
  {
    nentries = nevt;
  }
  for (int ientry=0; ientry<nentries; ientry++)
  {
    tree->GetEntry(ientry);

    if (ientry<4)
    {
      // print charge from channels 0 and 127
      cout << f_evt << "\t" << f_tt[0] << "\t" << f_tq[0] << "\t" << f_q[0] << endl;
      cout << "\t" << f_tt[127] << "\t" << f_tq[127] << "\t" << f_q[127] << endl;
    }

    for (int ipmt=0; ipmt<NUM_PMT; ipmt++)
    {
      int ifeech = (ipmt/8)*16 + 8 + ipmt%8;
      float tt = f_tt[ipmt];
      float tq = f_tq[ipmt];

      if ( pass>0 )
      {
        if ( !isnan(mcal->get_tt0(ipmt)) )
        {
          tt -= mcal->get_tt0(ipmt);
        }
        if ( !isnan(mcal->get_tq0(ipmt)) )
        {
          tq -= mcal->get_tq0(ipmt);
        }
      }

      h_tt[ipmt]->Fill( tt );
      h2_tt->Fill( tt, ipmt );

      h_tq[ipmt]->Fill( tq );
      h2_tq->Fill( tq, ipmt );

      int arm = ipmt/64;

      if ( pass==0 && tt>0. && tt<26. )
      {
        h_q[ipmt]->Fill( f_q[ipmt] );
      }
      else if ( pass==1 && fabs(tt)<26. )
      {
        h_q[ipmt]->Fill( f_q[ipmt] );
      }
      else if ( pass==2 && fabs(tq)<15 )
      {
        if ( arm==0 && f_bn[0]<30 )
        {
          h_q[ipmt]->Fill( f_q[ipmt] );
        }
        else if ( arm==1 && f_bn[1]<30 )
        {
          h_q[ipmt]->Fill( f_q[ipmt] );
        }
      }
    }
  }

  TCanvas *ac[100];
  int cvindex = 0;

  //== Calculate tq_t0 and tt_t0 on pass0, or apply tt_t0 if pass>0
  ac[cvindex] = new TCanvas("cal_tt","ch vs tt",425*1.5,550*1.5);
  h2_tt->Draw("colz");
  if ( pass==0 )
  {
    name = dir + "h2_tt.png";
  }
  else if ( pass>0 )
  {
    name = dir + "h2_ttcorr.png";
    h2_tt->GetXaxis()->SetRangeUser( 0.,20 );
    gPad->Modified();
    gPad->Update();
  }
  cout << name << endl;
  ac[cvindex]->Print( name );
  ++cvindex;

  // now tq_t0
  ac[cvindex] = new TCanvas("cal_tq","ch vs tq",425*1.5,550*1.5);
  h2_tq->Draw("colz");
  if ( pass==0 )
  {
    name = dir + "h2_tq.png";
  }
  else if ( pass>0 )
  {
    name = dir + "h2_tqcorr.png";
    h2_tq->GetXaxis()->SetRangeUser( -20,20 );
    gPad->Modified();
    gPad->Update();
  }
  cout << name << endl;
  ac[cvindex]->Print( name );
  ++cvindex;

  // Here we calculate tt_t0 and tq_t0, starting with tt_t0 first
  ac[cvindex] = new TCanvas("cal_tt_ch","tt",550*1.5,425*1.5);
  gPad->SetLogy(1);

  ofstream cal_tt_t0_file;
  if ( pass==0 ) 
  {
    name = dir + "bbc_tt_t0.calib";
    cal_tt_t0_file.open( name );
  }

  TF1 *gaussian = new TF1("gaussian","gaus",-25,25);
  gaussian->SetLineColor(2);
  double min_twindow = -25.;
  double max_twindow = 25.;

  for (int ipmt=0; ipmt<NUM_PMT && pass==0; ipmt++)
  {
    // only look in the middle
    if ( ipmt==0 || ipmt==64 )
    {
      // use wide range for 1st channel in each arm
      h_tt[ipmt]->SetAxisRange(0.,27.);
    }
    else
    {
      // for subsequent channels in an arm, use a window
      // determined by the 1st channel
      h_tt[ipmt]->SetAxisRange(min_twindow,max_twindow);
    }
    Double_t peak = h_tt[ipmt]->GetMaximum();
    int peakbin = h_tt[ipmt]->GetMaximumBin();
    Double_t mean = h_tt[ipmt]->GetBinCenter( peakbin );
    Double_t sigma = 1.0;
    gaussian->SetParameters( peak, mean, 5 );
    gaussian->SetRange( mean-3*sigma, mean+3*sigma );

    if ( ipmt==0 || ipmt==64 )
    {
      min_twindow = mean - 3*sigma;
      max_twindow = mean + 3*sigma;
    }

    h_tt[ipmt]->SetAxisRange(0.,27.);
    h_tt[ipmt]->Fit(gaussian,"R");
    //gPad->SetLogy(1);
    mean = gaussian->GetParameter(1);
    Double_t meanerr = gaussian->GetParError(1);
    sigma = gaussian->GetParameter(2);
    Double_t sigmaerr = gaussian->GetParError(2);

    if ( pass==0 )
    {
      cal_tt_t0_file << ipmt << "\t" << mean << "\t" << meanerr << "\t" << sigma << "\t" << sigmaerr << endl;
      name = dir + "h_tt"; name += ipmt; name += ".png";
    }
    else if ( pass==1 || pass==2 )
    {
      name = dir + "h_ttcorr"; name += ipmt; name += ".png";
    }
    cout << name << endl;
    ac[cvindex]->Print( name );
  }
  if ( pass==0 )
  {
    cal_tt_t0_file.close();
  }
  ++cvindex;

  // Now calculate tq_t0
  ac[cvindex] = new TCanvas("cal_tq_ch","tq",550*1.5,425*1.5);
  gPad->SetLogy(1);

  ofstream cal_tq_t0_file;
  if ( pass==0 ) 
  {
    name = dir + "bbc_tq_t0.calib";
    cal_tq_t0_file.open( name );
  }

  for (int ipmt=0; ipmt<NUM_PMT && pass==0; ipmt++)
  {
    // only look in the middle
    if ( ipmt==0 || ipmt==64 )
    {
      // use wide range for 1st channel in each arm
      h_tq[ipmt]->SetAxisRange(-25.,25.);
    }
    else
    {
      // for subsequent channels in an arm, use a window
      // determined by the 1st channel
      h_tq[ipmt]->SetAxisRange(min_twindow,max_twindow);
    }
    Double_t peak = h_tq[ipmt]->GetMaximum();
    int peakbin = h_tq[ipmt]->GetMaximumBin();
    Double_t mean = h_tq[ipmt]->GetBinCenter( peakbin );
    Double_t sigma = 1.0;
    gaussian->SetParameters( peak, mean, 5 );
    gaussian->SetRange( mean-3*sigma, mean+3*sigma );

    if ( ipmt==0 || ipmt==64 )
    {
      min_twindow = mean - 3*sigma;
      max_twindow = mean + 3*sigma;
    }

    h_tq[ipmt]->SetAxisRange(-25.,25.);
    h_tq[ipmt]->Fit(gaussian,"R");
    //gPad->SetLogy(1);
    mean = gaussian->GetParameter(1);
    Double_t meanerr = gaussian->GetParError(1);
    sigma = gaussian->GetParameter(2);
    Double_t sigmaerr = gaussian->GetParError(2);

    if ( pass==0 )
    {
      cal_tq_t0_file << ipmt << "\t" << mean << "\t" << meanerr << "\t" << sigma << "\t" << sigmaerr << endl;
      name = dir + "h_tq"; name += ipmt; name += ".png";
    }
    else if ( pass==1 || pass==2 )
    {
      name = dir + "h_tqcorr"; name += ipmt; name += ".png";
    }
    cout << name << endl;
    ac[cvindex]->Print( name );
  }
  if ( pass==0 )
  {
    cal_tq_t0_file.close();
  }
  ++cvindex;

  //== Draw the charge histograms
  ac[cvindex] = new TCanvas("cal_q","q",425*1.5,550*1.5);
  if ( pass==0 )
  {
    gPad->SetLogy(1);
  }

  for (int ipmt=0; ipmt<NUM_PMT && pass==0; ipmt++)
  {
      h_q[ipmt]->Draw();

      name = dir + "h_q"; name += ipmt; name += ".png";
      cout << name << endl;
      ac[cvindex]->Print( name );
  }
  ++cvindex;

  savefile->Write();
  //savefile->Close();

  if ( pass>0 )
  {
    recal_mbd_mip( savefname, pass, nevt );
  }
}

