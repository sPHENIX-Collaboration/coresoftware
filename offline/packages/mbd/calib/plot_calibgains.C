// Make plots of the values of the gains
#include "mbd/MbdCalib.h"
//#include "get_runstr.h"

R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libmbd_io.so)

void read_calibgains(const char *flist);
void plot_relwidth_onerun(const int irun = 0);

TFile *savefile;

const int MAX_RUNS = 1000;
int nruns = 0;

MbdCalib *bcal[MAX_RUNS];

Double_t bqmean[128][MAX_RUNS];
Double_t bqmeanerr[128][MAX_RUNS];
Double_t bqwidth[128][MAX_RUNS];
Double_t bqwidtherr[128][MAX_RUNS];
Double_t bqrelwidth[128][MAX_RUNS];
Double_t bqrelwidtherr[128][MAX_RUNS];
Double_t listofruns[MAX_RUNS];
Double_t runindex[MAX_RUNS];

TGraphErrors *gainvals[128];
TGraphErrors *g_relwidth[128];  // relative width
TH1 *h_relwidth{nullptr};

const int update_qfit = 0;    // whether to write update gain files
const int write_plots = 1;    // whether to save plots to png
  

// flist is a list of run-seq or run number directories
// where the calibrations are stored
void read_calibgains(const char *flist)
{
  ifstream inflist(flist);

  ifstream cal_mip_file;

  TString calrunseq;
  TString calfile;
  TString name;
  while ( inflist >> calrunseq )
  {
    calfile = "results/" + calrunseq + "/mbd_qfit.calib";
    //calfile = "results/" + calrunseq + "/bbc_qfit.calib";
    cout << calfile << endl;
    cal_mip_file.open( calfile );

    TString runtext = calrunseq;
    runtext.ReplaceAll("-0000","");
    listofruns[nruns] = runtext.Atof();
    cout << listofruns[nruns] << endl;

    runindex[nruns] = nruns;

    //bcal[nruns] = new MbdCalib();
    //bcal[nruns]->Download_Gains( calfile.Data() );

    int    temp_pmt;
    double integ;
    double best_peak;
    double width;
    double integerr;
    double best_peakerr;
    double widtherr;
    double chi2ndf;

    double corrected_peak;

    for (int ipmt=0; ipmt<128; ipmt++)
    {
      //float gain = bcal[nruns]->get_qgain(ipmt);
      //float gainerr = bcal[nruns]->get_qgainerr(ipmt);

      cal_mip_file >> temp_pmt >> integ >> best_peak >> width 
        >> integerr >> best_peakerr >> widtherr >> chi2ndf;

      /*
      cout << temp_pmt << "\t" << integ << "\t" << best_peak << "\t" << width 
        << "\t" << integerr << "\t" << best_peakerr << "\t" << widtherr << "\t" << chi2ndf << endl;
      */

      bqmean[temp_pmt][nruns] = best_peak;
      bqmeanerr[temp_pmt][nruns] = best_peakerr;
      bqwidth[temp_pmt][nruns] = fabs(width);
      bqwidtherr[temp_pmt][nruns] = widtherr;
      
      /*
      if ( listofruns[nruns] == 21520 )
      {
        cout << "XXX " << calfile << "\t" << ipmt << "\t" << best_peakerr << "\t"
          << best_peak << "\t" << bqmean[temp_pmt][nruns-1] << endl;
      }
      */

      // check for bad fit
      if ( integ < 0. || bqmean[temp_pmt][nruns]<0. )
      {
        cout << "BAD " << calfile << "\t" << ipmt << "\t" << integ << "\t" << best_peak << endl;
        bqmean[temp_pmt][nruns] = NAN;
      }

      if ( bqmeanerr[temp_pmt][nruns]>10. )
      {
        cout << "BADERR " << calfile << "\t" << ipmt << "\t" << best_peakerr << "\t"
          << best_peak << "\t" << bqmean[temp_pmt][nruns-1] << endl;
        bqmeanerr[temp_pmt][nruns] = NAN;
      }
    }

    cal_mip_file.close();

    nruns++;
  }
  inflist.close();
}

void plot_relwidth_onerun(const int irun = 0)
{
  if ( h_relwidth==nullptr )
  {
    h_relwidth = new TH1F("h_relwidth","relative width",100,0.1,0.36);
  }
  cout << "== Relative Widths ==" << endl;
  for (int ipmt=0; ipmt<128; ipmt++)
  {
    double relwidth = bqwidth[ipmt][irun]/bqmean[ipmt][irun];
    double frac_merr = bqmeanerr[ipmt][irun]/bqmean[ipmt][irun];
    double frac_werr = bqwidtherr[ipmt][irun]/bqwidth[ipmt][irun];
    double relwidtherr = frac_merr*frac_merr + frac_werr*frac_werr;
    relwidtherr = relwidth*sqrt(relwidtherr);

    cout << ipmt << "\t" << relwidth << "\t" << relwidtherr << endl;

    h_relwidth->Fill( relwidth );
  }
  h_relwidth->Draw();
}

void plot_calibgains(const char *flist = "runseq.list")
{
  // Read in all the calibrations from flist
  read_calibgains(flist);

  if ( write_plots )
  {
    savefile = new TFile("results/calibgains.root","RECREATE");
  }

  TString name;
  TString title;
  for (int ipmt=0; ipmt<128; ipmt++)
  {
    // do corrections
    for (int irun=0; irun<nruns; irun++)
    {
      if ( isnan(bqmean[ipmt][irun]) )
      {
        //if ( ipmt==51 && irun== )
        // search for next good val
        double nextgood = NAN;
        for (int nextrun=irun+1; nextrun<nruns; nextrun++)
        {
          if ( !isnan( bqmean[ipmt][nextrun] ) )
          {
            nextgood = bqmean[ipmt][nextrun];
            break;
          }
        }


        // search for prev good val
        double prevgood = NAN;
        for (int prevrun=irun-1; prevrun>=0; prevrun--)
        {
          if ( !isnan( bqmean[ipmt][prevrun] ) )
          {
            prevgood = bqmean[ipmt][prevrun];
            break;
          }
        }

        if ( !isnan(nextgood) && !isnan(prevgood) )
        {
          bqmean[ipmt][irun] = (nextgood+prevgood)/2.0;
        }
        else if ( !isnan(nextgood) )
        {
          bqmean[ipmt][irun] = nextgood;
        }
        else if ( !isnan(prevgood) )
        {
          bqmean[ipmt][irun] = prevgood;
        }
        else
        {
          cout << "ERROR, no good run to interpolate from" << endl;
        }
      }

      if ( isnan(bqmeanerr[ipmt][irun]) )
      {
        // search for next good val
        double nextgood = NAN;
        for (int nextrun=irun+1; nextrun<nruns; nextrun++)
        {
          if ( !isnan( bqmeanerr[ipmt][nextrun] ) )
          {
            nextgood = bqmeanerr[ipmt][nextrun];
            break;
          }
        }


        // search for prev good val
        double prevgood = NAN;
        for (int prevrun=irun-1; prevrun>=0; prevrun--)
        {
          if ( !isnan( bqmeanerr[ipmt][prevrun] ) )
          {
            prevgood = bqmeanerr[ipmt][prevrun];
            break;
          }
        }

        if ( irun==0 )
        {
          bqmeanerr[ipmt][irun] = nextgood;
        }
        else if ( irun==nruns-1 )
        {
          bqmeanerr[ipmt][irun] = prevgood;
        }
        else
        {
          bqmeanerr[ipmt][irun] = (nextgood+prevgood)/2.0;
        }
      }
    }

    name = "gainvals"; name += ipmt;
    title = "gain, ch"; title += ipmt;
    //gainvals[ipmt] = new TGraphErrors(nruns,runindex,bqmean[ipmt],0,bqmeanerr[ipmt]);
    gainvals[ipmt] = new TGraphErrors(nruns,listofruns,bqmean[ipmt],0,bqmeanerr[ipmt]);
    gainvals[ipmt]->SetName( name );
    gainvals[ipmt]->SetTitle( title );

    gainvals[ipmt]->Draw("ap");
    gPad->Modified();
    gPad->Update();

    if ( write_plots )
    {
      name += ".png";
      gPad->Print( name );
      gainvals[ipmt]->Write();
    }
  }

  // Get mean gain for each channel, and write out the values
  TString meangainfname = "results/"; meangainfname += flist;
  meangainfname.ReplaceAll(".list","_meangains.calib");
  ofstream meangainfile( meangainfname );
  TF1 *f_meangain[128] {nullptr};
  Double_t grp_mean[128];
  Double_t grp_meanerr[128];
  for (int ipmt=0; ipmt<128; ipmt++)
  {
    name = "f_meangain"; name += ipmt;
    f_meangain[ipmt] = new TF1(name,"pol0",0,1e9);
    gainvals[ipmt]->Fit(f_meangain[ipmt]);
    gPad->Update();
    gPad->Modified();
    grp_mean[ipmt] = f_meangain[ipmt]->GetParameter(0);
    grp_meanerr[ipmt] = f_meangain[ipmt]->GetParError(0);

    meangainfile << ipmt << "\t" << grp_mean[ipmt] << "\t" << grp_meanerr[ipmt] << endl;
    cout << ipmt << "\t" << grp_mean[ipmt] << "\t" << grp_meanerr[ipmt] << endl;
  }
  meangainfile.close();

  // generate new gains
  ifstream inflist;
  ifstream cal_mip_file;
  TString calrunseq;
  TString calfile;

  if ( update_qfit )
  {
    inflist.open(flist);
    ofstream newcal_mip_file;

    int    temp_pmt;
    double integ;
    double best_peak;
    double width;
    double integerr;
    double best_peakerr;
    double widtherr;
    double chi2ndf;

    int irun = 0;
    while ( inflist >> calrunseq )
    {
      calfile = "results/" + calrunseq + "/bbc_qfit.calib";
      cout << calfile << endl;
      cal_mip_file.open( calfile );

      // create new gain corr file
      calfile = "results/" + calrunseq + "/newgainvals_mbd_qfit.calib";
      newcal_mip_file.open( calfile );

      for (int ipmt=0; ipmt<128; ipmt++)
      {

        cal_mip_file >> temp_pmt >> integ >> best_peak >> width 
          >> integerr >> best_peakerr >> widtherr >> chi2ndf;

        newcal_mip_file << temp_pmt << "\t" << integ << "\t" << bqmean[ipmt][irun] << "\t" << width << "\t"
          << integerr << "\t" << bqmeanerr[ipmt][irun] << "\t" << widtherr << "\t"
          << chi2ndf << endl;
      }

      cal_mip_file.close();
      newcal_mip_file.close();
      irun++;
    }
  }

  if ( write_plots )
  {
    savefile->Write();
  }
}

