#include "MbdEvent.h"
#include "MbdCalib.h"
#include "MbdGeomV1.h"
#include "MbdOut.h"
#include "MbdPmtContainer.h"
#include "MbdPmtHit.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <sphenixnpc/CDBUtils.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>

#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>

#include <cmath>
#include <iomanip>
#include <iostream>

MbdEvent::MbdEvent()
  : _tres(0.05)
{
  // set default values

  int nsamples = 31;  /// NEED TO MAKE THIS FLEXIBLE
  recoConsts *rc = recoConsts::instance();
  if ( rc->FlagExist("MBD_TEMPLATEFIT") )
  {
    do_templatefit = rc->get_IntFlag("MBD_TEMPLATEFIT");
  }
  else
  {
    do_templatefit = 0;
  }

  for (int ifeech = 0; ifeech < MbdDefs::BBC_N_FEECH; ifeech++)
  {
    // std::cout << PHWHERE << "Creating _mbdsig " << ifeech << std::endl;
    _mbdsig.emplace_back(ifeech, nsamples);

    // Do evt-by-evt pedestal using sample range below
    //_mbdsig[ifeech].SetEventPed0Range(0,1);
    _mbdsig[ifeech].SetEventPed0PreSamp(6, 2);
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

  gaussian = nullptr;

  // BBCCALIB is used in offline to read in our calibrations
  const char *bbccaldir = getenv("BBCCALIB");
  if (bbccaldir)
  {
    // Online calibrations
    std::string gainfile = std::string(bbccaldir) + "/" + "bbc_mip.calib";
    Read_Charge_Calib( gainfile );

    std::string tq_t0_offsetfile = std::string(bbccaldir) + "/" + "bbc_tq_t0.calib";
    Read_TQ_T0_Offsets( tq_t0_offsetfile );

    std::string tq_clk_offsetfile = std::string(bbccaldir) + "/" + "bbc_tq_clk.calib";
    Read_TQ_CLK_Offsets( tq_clk_offsetfile );

    std::string tt_clk_offsetfile = std::string(bbccaldir) + "/" + "bbc_tt_clk.calib";
    Read_TT_CLK_Offsets( tt_clk_offsetfile );

    /*
    std::string mondata_fname = std::string(bbccaldir) + "/" + "BbcMonData.dat";
    std::ifstream mondatafile( mondata_fname );
    string label;
    mondatafile >> label >> bz_offset;
    std::cout << PHWHERE << label << "\t" << bz_offset << std::endl;
    mondatafile.close();
    */
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
  delete gaussian;
  delete _mbdgeom;
  delete _mbdcal;
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
  _mbdcal->Download_All();

  // Read in template if specified
  if ( do_templatefit )
  {
    for (int ifeech = 0; ifeech < MbdDefs::BBC_N_FEECH; ifeech++)
    {
      if ( _mbdgeom->get_type(ifeech) == 0 ) continue;
      //std::cout << PHWHERE << "Reading template " << ifeech << std::endl;
      //std::cout << "SIZES0 " << _mbdcal->get_shape(ifeech).size() << std::endl;
      _mbdsig[ifeech].SetTemplate( _mbdcal->get_shape(ifeech), _mbdcal->get_sherr(ifeech) );
      _mbdsig[ifeech].SetMinMaxFitTime( _mbdcal->get_sampmax(ifeech)-2-3, _mbdcal->get_sampmax(ifeech)-2+3 );
      //_mbdsig[ifeech].SetMinMaxFitTime( 0, 31 );
    }
  }

  return 0;
}

///
void MbdEvent::Clear()
{
  // Reset BBC/MBD raw data
  std::fill_n(m_pmttt, 128, 1e12);
  std::fill_n(m_pmttq, 128, 1e12);
  std::fill_n(m_pmtq, 128, 0.);

  // Reset BBC/MBD Arm Data
  for (int iarm = 0; iarm < 2; iarm++)
  {
    m_bbcn[iarm] = 0;
    m_bbcq[iarm] = 0.;
    m_bbct[iarm] = -9999.;
    m_bbcte[iarm] = -9999.;
    hevt_bbct[iarm]->Reset();
    hevt_bbct[1]->Reset();
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
      xmitclocks[ipkt] = static_cast<UShort_t>(p[ipkt]->iValue(0, "CLOCK"));

      femclocks[ipkt][0] = static_cast<UShort_t>(p[ipkt]->iValue(0, "FEMCLOCK"));
      femclocks[ipkt][1] = static_cast<UShort_t>(p[ipkt]->iValue(1, "FEMCLOCK"));

      for (int ich = 0; ich < NCHPERPKT; ich++)
      {
        int feech = ipkt * NCHPERPKT + ich;
        // std::cout << "feech " << feech << std::endl;
        for (int isamp = 0; isamp < MbdDefs::MAX_SAMPLES; isamp++)
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

        _mbdsig[feech].SetXY(m_samp[feech], m_adc[feech]);
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

  for (int ifeech = 0; ifeech < MbdDefs::BBC_N_FEECH; ifeech++)
  {
    int pmtch = _mbdgeom->get_pmt(ifeech);
    int type = _mbdgeom->get_type(ifeech);  // 0 = T-channel, 1 = Q-channel

    if (type == 0)
    {
      continue;
    }
    // Use dCFD method to get time for now in charge channels
    // Double_t threshold = 4.0*sig->GetPed0RMS();

    // std::cout << "getspline " << ifeech << std::endl;
    if ( do_templatefit )
    {
      _mbdsig[ifeech].FitTemplate();
 
      m_pmttq[pmtch] = _mbdsig[ifeech].GetTime();
      m_ampl[ifeech] = _mbdsig[ifeech].GetAmpl();
    }
    else
    {
      _mbdsig[ifeech].GetSplineAmpl();
      Double_t threshold = 0.5;
      m_pmttq[pmtch] = _mbdsig[ifeech].dCFD(threshold);
      m_ampl[ifeech] = _mbdsig[ifeech].GetAmpl();
    }

    if ( m_ampl[ifeech] < _mbdcal->get_qgain(pmtch)*0.25 )
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
    }

    m_pmtq[pmtch] = m_ampl[ifeech] / _mbdcal->get_qgain(pmtch);

    if (m_pmtq[pmtch] < 0.25)
    {
      m_pmtq[pmtch] = 0.;
      m_pmttq[pmtch] = std::numeric_limits<Float_t>::quiet_NaN();
    }

    // set tt to tq for now
    m_pmttt[pmtch] = m_pmttq[pmtch];

    /*
    if ( m_evt<3 && ifeech==255 && m_ampl[ifeech] )
    {
      std::cout << "dcfdcalc " << m_evt << "\t" << ifeech << "\t" << m_pmttq[pmtch] << "\t" << m_ampl[ifeech] << std::endl;
    }
    */
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
  return m_evt;
}

///
int MbdEvent::Calculate(MbdPmtContainer *bbcpmts, MbdOut *bbcout)
{
  if (_verbose >= 10)
  {
    std::cout << "In MbdEvent::Calculate()" << std::endl;
  }
  Clear();
  if (bbcout != nullptr)
  {
    bbcout->Reset();
  }

  if (!gaussian)
  {
    gaussian = new TF1("gaussian", "gaus", 0, 20);
    gaussian->FixParameter(2, _tres);  // set sigma to timing resolution
  }

  std::vector<float> hit_times[2];  // times of the hits in each [arm]

  // calculate bbc global variables
  if (_verbose >= 10)
  {
    std::cout << "Hit PMT info " << std::endl;
  }
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
    // std::cout << "earliest" << iarm << "\t" << earliest << std::endl;

    gaussian->SetParameter(0, 5);
    // gaussian->SetParameter(1, earliest);
    // gaussian->SetRange(6, earliest + 5 * 0.05);
    gaussian->SetParameter(1, hevt_bbct[iarm]->GetMean());
    gaussian->SetParameter(2, hevt_bbct[iarm]->GetRMS());
    gaussian->SetRange(hevt_bbct[iarm]->GetMean() - 5, hevt_bbct[iarm]->GetMean() + 5);

    if (_verbose)
    {
      if (ac == nullptr)
      {
        ac = new TCanvas("ac", "ac", 550 * 1.5, 425 * 1.5);
        ac->Divide(2, 1);
      }
      ac->cd(iarm + 1);
    }

    hevt_bbct[iarm]->Fit(gaussian, "BNQLR");
    if (_verbose)
    {
      hevt_bbct[iarm]->Draw();
    }

    // m_bbct[iarm] = m_bbct[iarm] / m_bbcn[iarm];
    m_bbct[iarm] = gaussian->GetParameter(1);
    m_bbcte[iarm] = earliest;

    //_bbcout->set_arm(iarm, m_bbcn[iarm], m_bbcq[iarm], m_bbct[iarm]);
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

    if (_verbose > 10)
    {
      std::cout << "bbcz " << m_bbcz << std::endl;
    }
  }

  // Fill rest of MbdOut
  if (bbcout != nullptr)
  {
    for (int iarm = 0; iarm < 2; iarm++)
    {
      bbcout->set_arm(iarm, get_bbcn(iarm), get_bbcq(iarm), get_bbct(iarm));
      // bbcout->set_clocks( m_evt, m_clk, m_femclk ); // only for V2
      if (_verbose > 10)
      {
        std::cout << get_bbcn(iarm) << "\t" << get_bbcq(iarm) << "\t" << get_bbct(iarm) << std::endl;
      }
    }

    if (get_bbcn(0) > 0 && get_bbcn(1) > 0)
    {
      bbcout->set_t0(get_bbct0(), get_bbct0err());
      bbcout->set_zvtx(get_bbcz(), get_bbczerr());
    }
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

