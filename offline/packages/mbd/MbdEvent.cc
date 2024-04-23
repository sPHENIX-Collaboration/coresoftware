#include "MbdEvent.h"
#include "MbdCalib.h"
#include "MbdGeomV1.h"
#include "MbdOut.h"
#include "MbdPmtContainer.h"
#include "MbdPmtHit.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>

#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TDirectory.h>

#include <cmath>
#include <iomanip>
#include <iostream>

MbdEvent::MbdEvent(const int cal_pass) :
  _tres(0.05),
  _calpass(cal_pass)
{
  // set default values

  _nsamples = MbdDefs::MAX_SAMPLES;  // Set to maximum initially, reset when we get a packet
  recoConsts *rc = recoConsts::instance();
  if (rc->FlagExist("MBD_TEMPLATEFIT"))
  {
    do_templatefit = rc->get_IntFlag("MBD_TEMPLATEFIT");
  }
  else
  {
    do_templatefit = 1;
  }

  for (int ifeech = 0; ifeech < MbdDefs::BBC_N_FEECH; ifeech++)
  {
    // std::cout << PHWHERE << "Creating _mbdsig " << ifeech << std::endl;
    _mbdsig.emplace_back(ifeech, _nsamples);
  }

  std::string name;
  std::string title;
  for (int iarm = 0; iarm < 2; iarm++)
  {
    //
    name = "hevt_bbct";
    name += std::to_string(iarm);
    title = "bbc times, arm ";
    title += std::to_string(iarm);
    hevt_bbct[iarm] = new TH1F(name.c_str(), title.c_str(), 2000, -50., 50.);
    hevt_bbct[iarm]->SetLineColor(4);
  }
  h2_tmax[0] = new TH2F("h2_ttmax", "time tmax vs ch", MbdDefs::MAX_SAMPLES, -0.5, MbdDefs::MAX_SAMPLES - 0.5, 128, 0, 128);
  h2_tmax[0]->SetXTitle("sample");
  h2_tmax[0]->SetYTitle("ch");
  h2_tmax[1] = new TH2F("h2_qtmax", "chg tmax vs ch", MbdDefs::MAX_SAMPLES, -0.5, MbdDefs::MAX_SAMPLES - 0.5, 128, 0, 128);
  h2_tmax[1]->SetXTitle("sample");
  h2_tmax[1]->SetYTitle("ch");

  for (float &iboard : TRIG_SAMP)
  {
    iboard = -1;
  }

  // Maybe this is used in the online monitor???
  // BBCCALIB is used in offline to read in our calibrations
  const char *bbccaldir = getenv("BBCCALIB");
  if (bbccaldir)
  {
    // Online calibrations
    std::string gainfile = std::string(bbccaldir) + "/" + "bbc_mip.calib";
    Read_Charge_Calib(gainfile);

    std::string tq_t0_offsetfile = std::string(bbccaldir) + "/" + "bbc_tq_t0.calib";
    Read_TQ_T0_Offsets(tq_t0_offsetfile);

    std::string tq_clk_offsetfile = std::string(bbccaldir) + "/" + "bbc_tq_clk.calib";
    Read_TQ_CLK_Offsets(tq_clk_offsetfile);

    std::string tt_clk_offsetfile = std::string(bbccaldir) + "/" + "bbc_tt_clk.calib";
    Read_TT_CLK_Offsets(tt_clk_offsetfile);

    /*
    std::string mondata_fname = std::string(bbccaldir) + "/" + "BbcMonData.dat";
    std::ifstream mondatafile( mondata_fname );
    string label;
    mondatafile >> label >> bz_offset;
    std::cout << PHWHERE << label << "\t" << bz_offset << std::endl;
    mondatafile.close();
    */
  }

  // Debug stuff
  _debugintt = 0;
  if (_debugintt )
  {
    ReadSyncFile();
  }

  Clear();
}

///
MbdEvent::~MbdEvent()
{
  for (auto &iarm : hevt_bbct)
  {
    delete iarm;
  }

  delete h2_tmax[0];
  delete h2_tmax[1];
  delete ac;
  delete gausfit[0];
  delete gausfit[1];
  delete _mbdgeom;
  delete _mbdcal;
  delete _syncttree;
}

int MbdEvent::InitRun()
{
  h2_tmax[0]->Reset();
  h2_tmax[1]->Reset();

  Clear();

  recoConsts *rc = recoConsts::instance();
  _runnum = rc->get_IntFlag("RUNNUMBER");
  if (_verbose)
  {
    std::cout << PHWHERE << "RUNNUMBER " << _runnum << std::endl;
  }

  if (_mbdgeom == nullptr)
  {
    _mbdgeom = new MbdGeomV1();
  }

  // Always reload calibrations on InitRun()
  if (_mbdcal != nullptr)
  {
    delete _mbdcal;
  }
  _mbdcal = new MbdCalib();
  std::cout << "SIMFLAG IS " << _simflag << std::endl;
  if (!_simflag)
  {
    _mbdcal->Download_All();
  }

  // Init parameters of the signal processing
  for (int ifeech = 0; ifeech < MbdDefs::BBC_N_FEECH; ifeech++)
  {
    // Do evt-by-evt pedestal using sample range below
    //_mbdsig[ifeech].SetEventPed0Range(0,1);
    //_mbdsig[ifeech].SetEventPed0PreSamp(6, 2, _mbdcal->get_sampmax(ifeech));
    if ( _calpass!=1 )
    {
      _mbdsig[ifeech].SetEventPed0PreSamp(5, 4, _mbdcal->get_sampmax(ifeech));
    }
    else
    {
      _mbdsig[ifeech].SetEventPed0Range(0,1);
    }

    // Read in template if specified
    if ( do_templatefit && _mbdgeom->get_type(ifeech)==1 )
    {
      // std::cout << PHWHERE << "Reading template " << ifeech << std::endl;
      // std::cout << "SIZES0 " << _mbdcal->get_shape(ifeech).size() << std::endl;
      //  Should set template size automatically here
      _mbdsig[ifeech].SetTemplate(_mbdcal->get_shape(ifeech), _mbdcal->get_sherr(ifeech));
      _mbdsig[ifeech].SetMinMaxFitTime(_mbdcal->get_sampmax(ifeech) - 2 - 3, _mbdcal->get_sampmax(ifeech) - 2 + 3);
      //_mbdsig[ifeech].SetMinMaxFitTime( 0, 31 );
    }
  }

  if ( _calpass > 0 )
  {
    _caldir = "results/"; _caldir += _runnum; _caldir += "/";
    TString cmd = "mkdir -p " + _caldir;
    gSystem->Exec( cmd );
    std::cout << "OUTPUT CALDIR = " << _caldir << std::endl;
  }

  if ( _calpass == 1 )
  {
    std::cout << "MBD Cal Pass 1" << std::endl;
    TDirectory *orig_dir = gDirectory;

    TString savefname = _caldir; savefname += "calpass1.root";
    std::cout << "Saving calpass 1 results to " << savefname << std::endl;
    _smax_tfile = std::make_unique<TFile>(savefname,"RECREATE");

    std::string name;

    for (int ich=0; ich<MbdDefs::MBD_N_FEECH; ich++)
    {
      name = "h_tmax"; name += std::to_string(ich);
      h_tmax[ich] = new TH1F(name.c_str(),name.c_str(),_nsamples,-0.5,_nsamples-0.5);
      h_tmax[ich]->SetXTitle("sample");
      h_tmax[ich]->SetYTitle("ch");
    }
    h2_tmax[0] = new TH2F("h2_ttmax","time tmax vs ch",_nsamples,-0.5,_nsamples-0.5,128,0,128);
    h2_tmax[1] = new TH2F("h2_qtmax","chg tmax vs ch",_nsamples,-0.5,_nsamples-0.5,128,0,128);
    h2_wave[0] = new TH2F("h2_twave","time adc vs ch",_nsamples,-0.5,_nsamples-0.5,128,0,128);
    h2_wave[1] = new TH2F("h2_qwave","chg adc vs ch",_nsamples,-0.5,_nsamples-0.5,128,0,128);

    h2_trange_raw = new TH2F("h2_trange_raw","tadc (raw) at trig samp vs ch",1600,0,16384,128,0,128);
    h2_trange = new TH2F("h2_trange","tadc at trig samp vs ch",1638,-100,16280,128,0,128);

    for (int itype=0; itype<2; itype++)
    {
      h2_tmax[itype]->SetXTitle("sample");
      h2_tmax[itype]->SetYTitle("ch");

      h2_wave[itype]->SetXTitle("sample");
      h2_wave[itype]->SetYTitle("ch");
    }

    orig_dir->cd();
  }
  else if ( _calpass == 2 )
  {
    // zero out the tt_t0, tq_t0, and gains to produce uncalibrated time and charge
    std::cout << "MBD Cal Pass 2" << std::endl;
    _mbdcal->Reset_TTT0();
    _mbdcal->Reset_TQT0();
    _mbdcal->Reset_Gains();
  }

  return 0;
}

int MbdEvent::End()
{
  if ( _calpass == 1 )
  {
    CalcSampMaxCalib();

    std::string fname = _caldir.Data(); fname += "mbd_sampmax.calib";
    _mbdcal->Write_SampMax( fname );

    _smax_tfile->Write();
  }

  return 1;
}

///
void MbdEvent::Clear()
{
  // Reset BBC/MBD raw data
  std::fill_n(m_pmttt, 128, std::numeric_limits<Float_t>::quiet_NaN());
  std::fill_n(m_pmttq, 128, std::numeric_limits<Float_t>::quiet_NaN());
  std::fill_n(m_pmtq, 128, 0.);

  // Reset BBC/MBD Arm Data
  for (int iarm = 0; iarm < 2; iarm++)
  {
    m_bbcn[iarm] = 0;
    m_bbcq[iarm] = 0.;
    m_bbct[iarm] = std::numeric_limits<Float_t>::quiet_NaN();
    m_bbcte[iarm] = std::numeric_limits<Float_t>::quiet_NaN();
    m_bbctl[iarm] = std::numeric_limits<Float_t>::quiet_NaN();
    hevt_bbct[iarm]->Reset();
    hevt_bbct[iarm]->GetXaxis()->SetRangeUser(-50, 50);
  }

  // Reset end product to prepare next event
  m_bbcz = std::numeric_limits<Float_t>::quiet_NaN();
  m_bbczerr = std::numeric_limits<Float_t>::quiet_NaN();
  m_bbct0 = std::numeric_limits<Float_t>::quiet_NaN();
  m_bbct0err = std::numeric_limits<Float_t>::quiet_NaN();
}

int MbdEvent::SetRawData(Event *event, MbdPmtContainer *bbcpmts)
{
  // First check if there is any event (ie, reading from PRDF)
  if (event == nullptr || bbcpmts == nullptr)
  {
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  int evt_type = event->getEvtType();
  if (evt_type != DATAEVENT)
  {
    std::cout << PHWHERE << "MbdEvent: Event type is not DATAEVENT, skipping" << std::endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  m_evt = event->getEvtSequence();
  UShort_t xmitclocks[2];    // [ipkt]
  UShort_t femclocks[2][2];  // [ipkt][iadc]

  // Get the relevant packets from the Event object and transfer the
  // data to the subsystem-specific table.

  // int flag_err = 0;
  for (int ipkt = 0; ipkt < 2; ipkt++)
  {
    int pktid = 1001 + ipkt;  // packet id
    p[ipkt] = event->getPacket(pktid);

    if (Verbosity() > 0)
    {
      static int counter = 0;
      if (counter < 4)
      {
        std::cout << "Found packet " << pktid << "\t" << p[ipkt] << std::endl;
        counter++;
      }
    }
    if (p[ipkt])
    {
      _nsamples = p[ipkt]->iValue(0, "SAMPLES");
      {
        static int counter = 0;
        if ( counter<2 )
        {
          std::cout << "NSAMPLES = " << _nsamples << std::endl;
        }
        counter++;
      }

      xmitclocks[ipkt] = static_cast<UShort_t>(p[ipkt]->iValue(0, "CLOCK"));

      femclocks[ipkt][0] = static_cast<UShort_t>(p[ipkt]->iValue(0, "FEMCLOCK"));
      femclocks[ipkt][1] = static_cast<UShort_t>(p[ipkt]->iValue(1, "FEMCLOCK"));

      for (int ich = 0; ich < NCHPERPKT; ich++)
      {
        int feech = ipkt * NCHPERPKT + ich;
        // std::cout << "feech " << feech << std::endl;
        for (int isamp = 0; isamp < _nsamples; isamp++)
        {
          m_adc[feech][isamp] = p[ipkt]->iValue(isamp, ich);
          m_samp[feech][isamp] = isamp;

          /*
          if ( m_adc[feech][isamp] <= 100 )
          {
            //flag_err = 1;
            //std::cout << "BAD " << m_evt << "\t" << feech << "\t" << m_samp[feech][isamp]
            //    << "\t" << m_adc[feech][isamp] << std::endl;
          }
          */
        }

        _mbdsig[feech].SetNSamples( _nsamples );
        _mbdsig[feech].SetXY(m_samp[feech], m_adc[feech]);
        //_mbdsig[feech].Print();
      }

      delete p[ipkt];
      p[ipkt] = nullptr;
    }
    else
    {
      // flag_err = 1;
      std::cout << PHWHERE << " ERROR, evt " << m_evt << " Missing Packet " << pktid << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  // Do a quick sanity check that all fem counters agree
  if (xmitclocks[0] != xmitclocks[1])
  {
    std::cout << __FILE__ << ":" << __LINE__ << " ERROR, xmitclocks don't agree" << std::endl;
  }
  for (auto &femclock : femclocks)
  {
    for (unsigned short iadc : femclock)
    {
      if (iadc != femclocks[0][0])
      {
        std::cout << __FILE__ << ":" << __LINE__ << " ERROR, femclocks don't agree" << std::endl;
      }
    }
  }

  // Store the clock info. We use just the first one, and assume all are consistent.
  m_clk = xmitclocks[0];
  m_femclk = femclocks[0][0];

  // We get SAMPMAX on this pass
  if ( _calpass == 1 )
  {
    FillSampMaxCalib();
    m_evt++;
    return -1001; // stop processing event (negative return values end event processing)
  }

  std::array<Double_t,MbdDefs::MBD_N_FEECH> tdc{0.};
  tdc.fill( 0. );

  for (int ifeech = 0; ifeech < MbdDefs::BBC_N_FEECH; ifeech++)
  {
    int pmtch = _mbdgeom->get_pmt(ifeech);
    int type = _mbdgeom->get_type(ifeech);  // 0 = T-channel, 1 = Q-channel

    // time channel
    if (type == 0)
    {
      tdc[pmtch] = _mbdsig[ifeech].MBDTDC(_mbdcal->get_sampmax(ifeech));

      if ( tdc[pmtch] < 40. || isnan(tdc[pmtch]) || fabs(_mbdcal->get_tt0(pmtch))>100. )
      {
        m_pmttt[pmtch] = std::numeric_limits<Float_t>::quiet_NaN();  // no hit
      }
      else
      {
        m_pmttt[pmtch] = _mbdcal->get_tcorr(ifeech,tdc[pmtch]);

        // at calpass 2, we use tcorr (uncal_mbd pass). make sure tt_t0 = 0.
        m_pmttt[pmtch] -= _mbdcal->get_tt0(pmtch);
      }

    }
    //else if ( type == 1 && !isnan(m_pmttt[pmtch]) ) // process charge channels which have good time hit
    else if ( type == 1 ) // process charge channels which have good time hit
    {

      // Use dCFD method to seed time in charge channels (or as primary if not fitting template)
      // std::cout << "getspline " << ifeech << std::endl;
      _mbdsig[ifeech].GetSplineAmpl();
      Double_t threshold = 0.5;
      m_pmttq[pmtch] = _mbdsig[ifeech].dCFD(threshold);
      m_ampl[ifeech] = _mbdsig[ifeech].GetAmpl(); // in adc units
      if (do_templatefit)
      {
        //std::cout << "fittemplate" << std::endl;
        _mbdsig[ifeech].FitTemplate( _mbdcal->get_sampmax(ifeech) );

        if ( _verbose )
        {
          std::cout << "tt " << ifeech << " " << pmtch << " " << m_pmttt[pmtch] << std::endl;
        }
        m_pmttq[pmtch] = _mbdsig[ifeech].GetTime(); // in units of sample number
        m_ampl[ifeech] = _mbdsig[ifeech].GetAmpl(); // in units of adc
      }

      // calpass 2, uncal_mbd. template fit. make sure qgain = 1, tq_t0 = 0
 
      if (m_ampl[ifeech] < _mbdcal->get_qgain(pmtch) * 0.25)
      {
        // m_t0[ifeech] = -9999.;
        m_pmttq[pmtch] = std::numeric_limits<Float_t>::quiet_NaN();
      }
      else
      {
        // if ( m_pmttq[pmtch]<-50. && ifeech==255 ) std::cout << "hit_times " << ifeech << "\t" << m_pmttq[pmtch] << std::endl;
        // if ( arm==1 ) std::cout << "hit_times " << ifeech << "\t" << setw(10) << m_pmttq[pmtch] << "\t" << board << "\t" << TRIG_SAMP[board] << std::endl;
        m_pmttq[pmtch] -= (_mbdcal->get_sampmax(ifeech) - 2);
        m_pmttq[pmtch] *= 17.7623;  // convert from sample to ns (1 sample = 1/56.299 MHz)
        m_pmttq[pmtch] = m_pmttq[pmtch] - _mbdcal->get_tq0(pmtch);

        // if tt is bad, use tq
        if ( fabs(_mbdcal->get_tt0(pmtch))>100. )
        {
          m_pmttt[pmtch] = m_pmttq[pmtch];
        }
        else
        {
          // we have a good tt ch. correct for slew if there is a hit
          //if ( ifeech==0 ) std::cout << "applying scorr" << std::endl;
          if ( !isnan(m_pmttt[pmtch]) )
          {
            m_pmttt[pmtch] -= _mbdcal->get_scorr(ifeech-8,m_ampl[ifeech]);
          }
        }
      }

      m_pmtq[pmtch] = m_ampl[ifeech] / _mbdcal->get_qgain(pmtch);

      if (m_pmtq[pmtch] < 0.25)
      {
        m_pmtq[pmtch] = 0.;
        m_pmttq[pmtch] = std::numeric_limits<Float_t>::quiet_NaN();
      }

      /*
      if ( m_evt<3 && ifeech==255 && m_ampl[ifeech] )
      {
        std::cout << "dcfdcalc " << m_evt << "\t" << ifeech << "\t" << m_pmttq[pmtch] << "\t" << m_ampl[ifeech] << std::endl;
      }
      */
    }
    else
    {
      m_pmtq[pmtch] = 0.;
      m_pmttq[pmtch] = std::numeric_limits<Float_t>::quiet_NaN();
    }
  }


  // bbcpmts->Reset();
  // std::cout << "q10 " << bbcpmts->get_tower_at_channel(10)->get_q() << std::endl;

  // Copy to output
  for (int ipmt = 0; ipmt < MbdDefs::BBC_N_PMT; ipmt++)
  {
    bbcpmts->get_pmt(ipmt)->set_pmt(ipmt, m_pmtq[ipmt], m_pmttt[ipmt], m_pmttq[ipmt]);
  }
  bbcpmts->set_npmt(MbdDefs::BBC_N_PMT);

  m_evt++;

  // Have uncalibrated charge and time at this pass
  if ( _calpass == 2 )
  {
    return -1002;
  }

  return m_evt;
}

///
int MbdEvent::Calculate(MbdPmtContainer *bbcpmts, MbdOut *bbcout)
{
  if ( _debugintt )
  {
    _verbose = 100;
  }
  //_verbose = 100;
 
  if (_verbose >= 10)
  {
    std::cout << "In MbdEvent::Calculate() " << m_evt << std::endl;
  }
  Clear();
  if (bbcout != nullptr)
  {
    bbcout->Reset();
  }

  // Debug stuff
  if ( _debugintt && (bbevt[_syncevt] != (m_evt - 1)))
  {
    _verbose = 0;
    return 1;
  }

  if (gausfit[0] == nullptr)
  {
    TString name;
    for (int iarm = 0; iarm < 2; iarm++)
    {
      name = "gausfit";
      name += iarm;
      gausfit[iarm] = new TF1(name, "gaus", 0, 20);
      gausfit[iarm]->FixParameter(2, _tres);  // set sigma to timing resolution
      gausfit[iarm]->SetLineColor(2);
    }
  }

  std::vector<float> hit_times[2];  // times of the hits in each [arm]

  // calculate bbc global variables
  if (_verbose >= 10)
  {
    std::cout << "Hit PMT info " << std::endl;
  }

  int epmt[2]{-1, -1};  // pmt of earliest time
  // int lpmt[2] {-1,-1};        // pmt of latest time
  double tepmt[2]{1e9, 1e9};    // earliest time
  double tlpmt[2]{-1e9, -1e9};  // latest time

  for (int ipmt = 0; ipmt < MbdDefs::BBC_N_PMT; ipmt++)
  {
    MbdPmtHit *bbcpmt = bbcpmts->get_pmt(ipmt);
    int arm = ipmt / 64;

    float t_pmt = bbcpmt->get_time();  // hit time of pmt
    float q_pmt = bbcpmt->get_q();     // charge in pmt

    if (_verbose >= 10)
    {
      std::cout << ipmt << "\t" << t_pmt << std::endl;
    }

    if (fabs(t_pmt) < 25. && q_pmt > 0.)
    {
      hit_times[arm].push_back(t_pmt);
      hevt_bbct[arm]->Fill(t_pmt);

      m_bbcn[arm]++;
      m_bbcq[arm] += q_pmt;

      if (_verbose >= 10)
      {
        std::cout << ipmt << "\t" << t_pmt << "\t" << q_pmt << std::endl;

        if (t_pmt < tepmt[arm])
        {
          epmt[arm] = ipmt;
          tepmt[arm] = t_pmt;
        }
        if (t_pmt > tlpmt[arm])
        {
          // lpmt[arm] = ipmt;
          tlpmt[arm] = t_pmt;
        }
      }
    }
  }

  if (_verbose >= 10)
  {
    std::cout << "nhits " << m_bbcn[0] << "\t" << m_bbcn[1] << std::endl;
  }
  // std::cout << "bbcte " << m_bbcte[0] << "\t" << m_bbcte[1] << std::endl;

  for (int iarm = 0; iarm < 2; iarm++)
  {
    if (hit_times[iarm].empty())
    {
      // std::cout << "hit_times size == 0" << std::endl;
      continue;
    }

    // std::cout << "EARLIEST " << iarm << std::endl;
    // std::cout << "ERROR hit_times size == " << hit_times[iarm].size() << std::endl;

    std::sort(hit_times[iarm].begin(), hit_times[iarm].end());
    float earliest = hit_times[iarm].at(0);
    float latest = hit_times[iarm].back();
    // std::cout << "earliest" << iarm << "\t" << earliest << std::endl;

    gausfit[iarm]->SetParameter(0, 5);
    // gausfit[iarm]->SetParameter(1, earliest);
    // gausfit[iarm]->SetRange(6, earliest + 5 * 0.05);
    gausfit[iarm]->SetParameter(1, hevt_bbct[iarm]->GetMean());
    gausfit[iarm]->SetParameter(2, hevt_bbct[iarm]->GetRMS());
    gausfit[iarm]->SetRange(hevt_bbct[iarm]->GetMean() - 5, hevt_bbct[iarm]->GetMean() + 5);

    if (_verbose)
    {
      if (ac == nullptr)
      {
        ac = new TCanvas("ac", "ac", 550 * 1.5, 425 * 1.5);
        ac->Divide(2, 1);
      }
      ac->cd(iarm + 1);
    }

    hevt_bbct[iarm]->Fit(gausfit[iarm], "BNQLR");

    // m_bbct[iarm] = m_bbct[iarm] / m_bbcn[iarm];
    m_bbct[iarm] = gausfit[iarm]->GetParameter(1);
    m_bbcte[iarm] = earliest;
    m_bbctl[iarm] = latest;

    //_bbcout->set_arm(iarm, m_bbcn[iarm], m_bbcq[iarm], m_bbct[iarm]);

    // if ( _verbose && mybbz[_syncevt]< -40. )
    if (_verbose)
    {
      hevt_bbct[iarm]->GetXaxis()->SetRangeUser(tepmt[iarm] - 3., tlpmt[iarm] + 3.);
      // hevt_bbct[iarm]->GetXaxis()->SetRangeUser(-20,20);
      hevt_bbct[iarm]->Draw();
      gausfit[iarm]->Draw("same");
      gPad->Modified();
      gPad->Update();
      if (iarm == 1)
      {
        double zearly = (tepmt[0] - tepmt[1]) * MbdDefs::C / 2.0;
        double znew = (m_bbct[0] - m_bbct[1]) * MbdDefs::C / 2.0;

        if (_debugintt)
        {
          double intzdiff = intz[_syncevt] / 10. - mybbz[_syncevt];
          double intzediff = intz[_syncevt] / 10. - zearly;
          if (fabs(znew - mybbz[_syncevt]) > 0.1)
          {
            std::cout << "**ERR** " << znew << "\t" << mybbz[_syncevt] << std::endl;
          }
          std::string junk;
          std::cout << m_evt << "\t" << bbevt[_syncevt] << "\t" << m_bbct[0] << "\t" << m_bbct[1] << std::endl;
          std::cout << m_evt << " gmean " << gausfit[0]->GetParameter(1) << "\t" << gausfit[1]->GetParameter(1) << std::endl;
          std::cout << m_evt << " mean " << hevt_bbct[0]->GetMean(1) << "\t" << hevt_bbct[1]->GetMean(1) << std::endl;
          std::cout << m_evt << " gsigma " << gausfit[0]->GetParameter(2) << "\t" << gausfit[1]->GetParameter(2) << std::endl;
          std::cout << m_evt << " rms " << hevt_bbct[0]->GetRMS() << "\t" << hevt_bbct[1]->GetRMS() << std::endl;
          std::cout << m_evt << " te ch " << epmt[0] << "\t" << epmt[1] << "\t" << tepmt[0] << "\t" << tepmt[1] << std::endl;
          std::cout << m_evt << " tetl " << m_bbcte[0] << "\t" << m_bbctl[0] << "\t" << m_bbcte[1] << "\t" << m_bbctl[1] << std::endl;
          std::cout << m_evt << " bz intz " << mybbz[_syncevt] << "\t" << intz[_syncevt] / 10. << "\t" << intzdiff << "\t" << intzdiff * 2.0 / MbdDefs::C << std::endl;
          std::cout << m_evt << " bze " << zearly << "\t" << intzediff << std::endl;
          std::cout << "? ";
          //std::cin >> junk;
        }
      }
    }
  }

  // Get Zvertex, T0
  if (m_bbcn[0] > 0 && m_bbcn[1] > 0)
  {
    // Now calculate zvtx, t0 from best times
    if (_verbose >= 10)
    {
      std::cout << "Evt " << m_evt << "\tt0\t" << m_bbct[0] << "\t" << m_bbct[1] << std::endl;
      std::cout << "bbcn " << m_bbcn[0] << "\t" << m_bbcn[1] << std::endl;
      std::cout << "bbcq " << m_bbcq[0] << "\t" << m_bbcq[1] << std::endl;
    }
    m_bbcz = (m_bbct[0] - m_bbct[1]) * TMath::C() * 1e-7 / 2.0;  // in cm
    m_bbct0 = (m_bbct[0] + m_bbct[1]) / 2.0;

    // correct z-vertex
    // m_bbcz += bz_offset;

    // hard code these for now
    // need study to determine muliplicity dependence
    m_bbczerr = 1.0;    // cm
    m_bbct0err = 0.05;  // ns

    /*
    // Use earliest time
    //cout << "t0\t" << m_bbct[0] << "\t" << m_bbct[1] << std::endl;
    //cout << "te\t" << m_bbcte[0] << "\t" << m_bbcte[1] << std::endl;
    m_bbcz = (m_bbcte[0] - m_bbcte[1]) * TMath::C() * 1e-7 / 2.0; // in cm
    m_bbct0 = (m_bbcte[0] + m_bbcte[1]) / 2.0;
    */

    // if (_verbose > 10)
    // if ( _verbose && mybbz[_syncevt]< -40. )
    if (_verbose)
    {
      std::cout << "bbcz " << m_bbcz << std::endl;
      std::string junk;
      std::cout << "? ";
      //std::cin >> junk;
    }
  }

  // Fill rest of MbdOut
  if (bbcout != nullptr)
  {
    for (int iarm = 0; iarm < 2; iarm++)
    {
      bbcout->set_arm(iarm, get_bbcn(iarm), get_bbcq(iarm), get_bbct(iarm));
      bbcout->set_clocks(m_evt, m_clk, m_femclk);  // only for V2
      if (_verbose > 10)
      {
        std::cout << get_bbcn(iarm) << "\t" << get_bbcq(iarm) << "\t" << get_bbct(iarm) << std::endl;
      }
    }

    if (get_bbcn(0) > 0 && get_bbcn(1) > 0)
    {
      bbcout->set_t0(get_bbct0(), get_bbct0err());
      bbcout->set_zvtx(get_bbcz(), get_bbczerr());
      
      if ( _debugintt )
      {
        bbcout->set_t0(intz[_syncevt]/10.);
      }

    }
  }

  if ( _debugintt )
  {
    _syncevt++;
    _verbose = 0;
  }

  return 1;
}

// This needs to be reconsidered for 2024 run, hopefully timing instability is fixed by then!
// Only used in online monitoring
int MbdEvent::DoQuickClockOffsetCalib()
{
  for (int ifeech = 0; ifeech < 256; ifeech++)
  {
    _mbdsig[ifeech].SetXY(m_samp[ifeech], m_adc[ifeech]);

    // determine the trig_samp board by board
    int tq = (ifeech / 8) % 2;  // 0 = T-channel, 1 = Q-channel
    int pmtch = (ifeech / 16) * 8 + ifeech % 8;

    double x, y;
    _mbdsig[ifeech].LocMax(x, y);
    h2_tmax[tq]->Fill(x, pmtch);
  }

  if (h2_tmax[1]->GetEntries() == 128 * 100)
  {
    TString name;
    TH1 *h_trigsamp[16]{};
    for (int iboard = 0; iboard < 16; iboard++)
    {
      name = "h_trigsamp";
      name += iboard;
      h_trigsamp[iboard] = h2_tmax[1]->ProjectionX(name, iboard * 8 + 1, (iboard + 1) * 8);
      int maxbin = h_trigsamp[iboard]->GetMaximumBin();
      TRIG_SAMP[iboard] = h_trigsamp[iboard]->GetBinCenter(maxbin);
      // std::cout << "iboard " << iboard << "\t" << iboard*8+1 << "\t" << (iboard+1)*8 << "\t" << h_trigsamp[iboard]->GetEntries() << std::endl;
      std::cout << "TRIG_SAMP" << iboard << "\t" << TRIG_SAMP[iboard] << std::endl;
    }
  }

  return 1;
}

// Store data for sampmax calibration (to correct ADC sample offsets by channel)
// This will replace DoQuickClockOffsetCalib() for 2024
// We might have to set this so it stops the sampmax calibration after N events
int MbdEvent::FillSampMaxCalib()
{
  for (int ifeech = 0; ifeech<MbdDefs::MBD_N_FEECH; ifeech++)
  {
    // determine the trig_samp board by board
    int type = _mbdgeom->get_type(ifeech);  // 0 = T-channel, 1 = Q-channel
    int pmtch = _mbdgeom->get_pmt(ifeech);
                                                                                  
    _mbdsig[ifeech].SetXY(m_samp[ifeech], m_adc[ifeech]);

    for (int isamp=0; isamp<_nsamples; isamp++)
    {
      // sanity check
      if ( m_samp[ifeech][isamp] != isamp )
      {
        std::cerr << PHWHERE << ", ch" << ifeech << ", msamp != isamp, " << m_samp[ifeech][isamp] << " " << isamp << std::endl;
      }
      h2_wave[type]->Fill( isamp , pmtch, m_adc[ifeech][isamp] );
    }

    double maxsamp, maxadc;
    _mbdsig[ifeech].LocMax( maxsamp, maxadc );

    h_tmax[ifeech]->Fill( maxsamp );
    h2_tmax[type]->Fill( maxsamp, pmtch );
    //std::cout << "fillint h2_tmax " << pmtch << "\t" << maxsamp << std::endl;
    //_mbdsig[ifeech].Print();

    if ( type==0 && h2_tmax[0]->GetEntries() > (128*1000) )
    {
      int samp_max = _mbdcal->get_sampmax( ifeech );

      h2_trange_raw->Fill( m_adc[ifeech][samp_max] , pmtch );

      TGraphErrors *gsubpulse = _mbdsig[ifeech].GetGraph();
      Double_t *y = gsubpulse->GetY();
      Double_t tdc = y[samp_max];
      h2_trange->Fill( tdc , pmtch );
    }

  }

  // At 1000 events, get the tmax for the time channels
  // so we can fill the h2_trange histograms
  if ( h2_tmax[0]->GetEntries() == (128*1000) )
  {
    CalcSampMaxCalib();
  }

  return 1;
}

int MbdEvent::CalcSampMaxCalib()
{
  TDirectory *orig_dir = gDirectory;
  _smax_tfile->cd();

  // sampmax for each board, for time and ch channels
  int feech = 0;
  for (int iboard=0; iboard<16; iboard++)
  {
    int min_ybin = iboard*8 + 1;
    int max_ybin = iboard*8 + 8;

    // sampmax for time channels
    TString name = "t_sampmax_bd"; name += iboard;
    TH1 *h_projx = h2_tmax[0]->ProjectionX(name,min_ybin,max_ybin);
    int maxbin = h_projx->GetMaximumBin();
    int samp_max = h_projx->GetBinCenter( maxbin );
    for (int ich=0; ich<8; ich++)
    {
      _mbdcal->set_sampmax( feech, samp_max );
      feech++;
    }
    delete h_projx;

    // sampmax for charge channels
    name = "t_sampmax_bd"; name += iboard;
    h_projx = h2_tmax[1]->ProjectionX(name,min_ybin,max_ybin);
    maxbin = h_projx->GetMaximumBin();
    samp_max = h_projx->GetBinCenter( maxbin );
    for (int ich=0; ich<8; ich++)
    {
      _mbdcal->set_sampmax( feech, samp_max );
      feech++;
    }
    delete h_projx;
  }

  orig_dir->cd();

  return 1;
}

int MbdEvent::Read_Charge_Calib(const std::string &gainfname)
{
  std::ifstream gainfile(gainfname);

  std::cout << "Reading gains from " << gainfname << std::endl;
  int ch;
  float integ, integerr;
  float peak, peakerr;
  float width, widtherr;
  float chi2ndf;
  while (gainfile >> ch >> integ >> peak >> width >> integerr >> peakerr >> widtherr >> chi2ndf)
  {
    gaincorr[ch] = 1.0 / peak;

    // std::cout << ch << "\t" << peak << std::endl;
  }

  gainfile.close();

  return 1;
}

// Read in tq t0 offset calibrations
int MbdEvent::Read_TQ_T0_Offsets(const std::string &t0cal_fname)
{
  std::ifstream tcalibfile(t0cal_fname);

  std::cout << "Reading tq_t0 offset calibrations from " << t0cal_fname << std::endl;

  int pmtnum;
  float meanerr;
  float sigma;
  float sigmaerr;
  for (int ipmt = 0; ipmt < MbdDefs::BBC_N_PMT; ipmt++)
  {
    tcalibfile >> pmtnum >> tq_t0_offsets[ipmt] >> meanerr >> sigma >> sigmaerr;
    if (pmtnum != ipmt)
    {
      std::cout << "ERROR, pmtnum != ipmt, " << pmtnum << "\t" << ipmt << std::endl;
    }
  }

  tcalibfile.close();

  return 1;
}

// Read in tq clk offset calibrations
int MbdEvent::Read_TQ_CLK_Offsets(const std::string &t0cal_fname)
{
  std::ifstream tcalibfile(t0cal_fname);

  std::cout << "Reading tq_clk offset calibrations from " << t0cal_fname << std::endl;

  int pmtnum;
  for (int ipmt = 0; ipmt < MbdDefs::BBC_N_PMT; ipmt++)
  {
    tcalibfile >> pmtnum >> tq_clk_offsets[ipmt];
    if (pmtnum != ipmt)
    {
      std::cout << "ERROR, pmtnum != ipmt, " << pmtnum << "\t" << ipmt << std::endl;
    }
  }

  tcalibfile.close();

  return 1;
}

// Read in tt clk offset calibrations
int MbdEvent::Read_TT_CLK_Offsets(const std::string &t0cal_fname)
{
  std::ifstream tcalibfile(t0cal_fname);

  std::cout << "Reading tq_clk offset calibrations from " << t0cal_fname << std::endl;

  int pmtnum;
  for (int ipmt = 0; ipmt < MbdDefs::BBC_N_PMT; ipmt++)
  {
    tcalibfile >> pmtnum >> tt_clk_offsets[ipmt];
    if (pmtnum != ipmt)
    {
      std::cout << "ERROR, pmtnum != ipmt, " << pmtnum << "\t" << ipmt << std::endl;
    }
  }

  tcalibfile.close();

  return 1;
}

void MbdEvent::ReadSyncFile(const char *fname)
{
  Int_t f_evt{0};
  UShort_t f_femclk{0};
  Float_t f_bz{0.};
  Long64_t bco_full{0};
  Double_t ES_zvtx{0.};
  Double_t mbd_bz{0.};

  _synctfile = std::make_unique<TFile>(fname, "READ");
  _syncttree = (TTree *) _synctfile->Get("t2");
  _syncttree->SetBranchAddress("evt", &f_evt);
  _syncttree->SetBranchAddress("femclk", &f_femclk);
  _syncttree->SetBranchAddress("bz", &f_bz);
  _syncttree->SetBranchAddress("bco_full", &bco_full);
  _syncttree->SetBranchAddress("ES_zvtx", &ES_zvtx);
  _syncttree->SetBranchAddress("mbd_bz", &mbd_bz);

  Stat_t nentries = _syncttree->GetEntries();
  for (int ientry = 0; ientry < nentries; ientry++)
  {
    _syncttree->GetEntry(ientry);

    bbevt.push_back(f_evt);
    bbclk.push_back(f_femclk);
    mybbz.push_back(f_bz);
    bco.push_back(bco_full);
    intz.push_back(ES_zvtx);
    bbz.push_back(mbd_bz);
  }

  std::cout << "Read in " << bbevt.size() << " INTT sync events" << std::endl;
}
