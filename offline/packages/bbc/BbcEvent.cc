#include "BbcEvent.h"
#include "BbcPmtInfoContainerV1.h"
#include "BbcOut.h"
#include "BbcGeomV1.h"
#include "BbcCalib.h"

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <phool/recoConsts.h>
#include <sphenixnpc/CDBUtils.h>

#include <TRandom.h>
#include <TString.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace BbcDefs;

BbcEvent::BbcEvent(void) :
  _verbose(0),
  _runnum(0),
  _eventnum(0),
  p{nullptr,nullptr},
  _tres(0.05),
  ac(nullptr)
{
  //set default values

  int nsamples = 31;  /// NEED TO MAKE THIS FLEXIBLE
  for (int ifeech=0; ifeech<BBC_N_FEECH; ifeech++)
  {
    //cout << "Creating bbcsig " << ifeech << endl;
    bbcsig.push_back( BbcSig(ifeech,nsamples) );
    
    // Do evt-by-evt pedestal using sample range below
    //bbcsig[ifeech].SetEventPed0Range(0,1);
    bbcsig[ifeech].SetEventPed0PreSamp(6,2);
  }

  TString name, title;
  for (int iarm = 0; iarm < 2; iarm++)
  {
    //
    name = "hevt_bbct";
    name += iarm;
    title = "bbc times, arm ";
    title += iarm;
    hevt_bbct[iarm] = new TH1F(name, title, 2000, -50., 50.);
    hevt_bbct[iarm]->SetLineColor(4);
  }
  h2_tmax[0] = new TH2F("h2_ttmax","time tmax vs ch",MAX_SAMPLES,-0.5,MAX_SAMPLES-0.5,128,0,128);
  h2_tmax[0]->SetXTitle("sample");
  h2_tmax[0]->SetYTitle("ch");
  h2_tmax[1] = new TH2F("h2_qtmax","chg tmax vs ch",MAX_SAMPLES,-0.5,MAX_SAMPLES-0.5,128,0,128);
  h2_tmax[1]->SetXTitle("sample");
  h2_tmax[1]->SetYTitle("ch");

  for (int iboard=0; iboard<16; iboard++)
  {
    TRIG_SAMP[iboard] = -1;
  }

  gaussian = nullptr;

  // BBCCALIB is used in offline to read in our calibrations 
  const char *bbccaldir = getenv("BBCCALIB");
  if ( bbccaldir )
  {
    // Online calibrations
    std::string gainfile = std::string(bbccaldir) + "/" + "bbc_mip.calib";
    Read_Charge_Calib( gainfile.c_str() );

    std::string tq_t0_offsetfile = std::string(bbccaldir) + "/" + "bbc_tq_t0.calib";
    Read_TQ_T0_Offsets( tq_t0_offsetfile.c_str() );

    std::string tq_clk_offsetfile = std::string(bbccaldir) + "/" + "bbc_tq_clk.calib";
    Read_TQ_CLK_Offsets( tq_clk_offsetfile.c_str() );

    std::string tt_clk_offsetfile = std::string(bbccaldir) + "/" + "bbc_tt_clk.calib";
    Read_TT_CLK_Offsets( tt_clk_offsetfile.c_str() );

    /*
    std::string mondata_fname = std::string(bbccaldir) + "/" + "BbcMonData.dat";
    ifstream mondatafile( mondata_fname );
    string label;
    mondatafile >> label >> bz_offset;
    std::cout << label << "\t" << bz_offset << std::endl;
    mondatafile.close();
    */
  }

  Clear();

}

///
BbcEvent::~BbcEvent()
{
  for (int iarm=0; iarm<2; iarm++)
  {
    delete hevt_bbct[iarm];
  }

  delete h2_tmax[0];
  delete h2_tmax[1];
  delete ac;
  delete gaussian;
  delete _bbcgeom;
  delete _bbccal;

}

int BbcEvent::InitRun()
{
  h2_tmax[0]->Reset();
  h2_tmax[1]->Reset();

  Clear();

  recoConsts *rc = recoConsts::instance();
  _runnum = rc->get_IntFlag("RUNNUMBER");
  if ( _verbose ) cout << "RUNNUMBER " << _runnum << endl;
  
  if ( _bbcgeom == nullptr )
  {
    _bbcgeom = new BbcGeomV1();
  }

  // Always reload calibrations on InitRun()
  if ( _bbccal != nullptr )
  {
    delete _bbccal;
  }
  _bbccal = new BbcCalib();
  _bbccal->Download_All();

  return 0;
}

///
void BbcEvent::Clear()
{
  // Reset BBC/MBD raw data
  std::fill_n(m_pmttt, 128, 1e12);
  std::fill_n(m_pmttq, 128, 1e12);
  std::fill_n(m_pmtq, 128, 0.);

  // Reset BBC/MBD Arm Data
  for ( int iarm = 0; iarm < 2; iarm++ )
  {
    m_bbcn[iarm] = 0;
    m_bbcq[iarm] = 0.;
    m_bbct[iarm] = -9999.;
    m_bbcte[iarm] = -9999.;
    hevt_bbct[iarm]->Reset();
    hevt_bbct[1]->Reset();
  }

  // Reset end product to prepare next event
  m_bbcz = NAN;
  m_bbczerr = NAN;
  m_bbct0 = NAN;
  m_bbct0err = NAN;
}


int BbcEvent::SetRawData(Event *event, BbcPmtInfoContainerV1 *bbcpmts)
{
  // First check if there is any event (ie, reading from PRDF)
  if ( event==0 || event==nullptr ) return -1;

  int evt_type = event->getEvtType();
  if ( evt_type != DATAEVENT )
  {
    cout << "BbcEvent: Event type is not DATAEVENT, skipping" << endl;
    return -2;
  }

  // Get the relevant packets from the Event object and transfer the
  // data to the subsystem-specific table.

  //int flag_err = 0;
  for (int ipkt=0; ipkt<2; ipkt++)
  {
    int pktid = 1001 + ipkt;    // packet id
    p[ipkt] = event->getPacket( pktid );

    static int counter = 0;
    if ( counter<4 )
    {
      cout << "Found packet " << pktid << "\t" << p[ipkt] << endl;
      counter++;
    }

    if ( p[ipkt] )
    {
      for (int ich=0; ich<NCHPERPKT; ich++)
      {
        int feech = ipkt*NCHPERPKT + ich;
        //cout << "feech " << feech << endl;
        for (int isamp=0; isamp<BbcDefs::MAX_SAMPLES; isamp++)
        {
          m_adc[feech][isamp] = p[ipkt]->iValue(isamp,ich);
          m_samp[feech][isamp] = isamp;

          /*
          if ( m_adc[feech][isamp] <= 100 )
          {
            //flag_err = 1;
            //cout << "BAD " << _eventnum << "\t" << feech << "\t" << m_samp[feech][isamp]
            //    << "\t" << m_adc[feech][isamp] << endl;
          }
          */
        }

        bbcsig[feech].SetXY(m_samp[feech],m_adc[feech]);
      }

      delete p[ipkt];
      p[ipkt] = nullptr;
    }
    else
    {
      //flag_err = 1;
      cout << "ERROR, evt " << _eventnum << " Missing Packet " << pktid << endl;
    }
  }

  for (int ifeech=0; ifeech<BBC_N_FEECH; ifeech++)
  {
    int pmtch = _bbcgeom->get_pmt(ifeech);
    int tq = _bbcgeom->get_type(ifeech);   // 0 = T-channel, 1 = Q-channel

    if ( tq == 0 ) continue;
    // Use dCFD method to get time for now in charge channels
    //Double_t threshold = 4.0*sig->GetPed0RMS();

    //cout << "getspline " << ifeech << endl;
    bbcsig[ifeech].GetSplineAmpl();
    Double_t threshold = 0.5;
    m_pmttq[pmtch] = bbcsig[ifeech].dCFD( threshold );
    m_ampl[ifeech] = bbcsig[ifeech].GetAmpl();

    if ( m_ampl[ifeech]<24 || ( fabs( m_pmttq[pmtch] ) >25 ) )
    {
      //m_t0[ifeech] = -9999.;
      m_pmttq[pmtch] = -9999.;
    }
    else
    {
      //if ( m_pmttq[pmtch]<-50. && ifeech==255 ) cout << "hit_times " << ifeech << "\t" << m_pmttq[pmtch] << endl;
      //if ( arm==1 ) cout << "hit_times " << ifeech << "\t" << setw(10) << m_pmttq[pmtch] << "\t" << board << "\t" << TRIG_SAMP[board] << endl;
      m_pmttq[pmtch] -= (_bbccal->get_sampmax(ifeech) - 2);
      m_pmttq[pmtch] *= 17.7623;               // convert from sample to ns (1 sample = 1/56.299 MHz)
      m_pmttq[pmtch] = m_pmttq[pmtch] - _bbccal->get_tq0(pmtch);
    }

    m_pmtq[pmtch] = m_ampl[ifeech] * _bbccal->get_qgain(pmtch);

    /*
    if ( _eventnum<3 && ifeech==255 && m_ampl[ifeech] )
    {
      cout << "dcfdcalc " << _eventnum << "\t" << ifeech << "\t" << m_pmttq[pmtch] << "\t" << m_ampl[ifeech] << endl;
    }
    */
  }

  //bbcpmts->Reset();
  //cout << "q10 " << bbcpmts->get_tower_at_channel(10)->get_q() << endl;

  // Copy to output
  for (int ipmt=0; ipmt<BBC_N_PMT; ipmt++)
  {
    bbcpmts->get_tower_at_channel(ipmt)->set_pmt(ipmt, m_pmtq[ipmt], m_pmttt[ipmt], m_pmttq[ipmt]);
  }

  _eventnum++;
  return _eventnum;
}

///
int BbcEvent::Calculate(BbcPmtInfoContainerV1 *bbcpmts, BbcOut *bbcout)
{
  if ( _verbose>=10 ) cout << "In BbcEvent::Calculate()" << endl;
  Clear();
  if ( bbcout!=0 ) bbcout->Reset();

  if ( ! gaussian )
  {
    gaussian = new TF1("gaussian", "gaus", 0, 20);
    gaussian->FixParameter(2, _tres);  // set sigma to timing resolution
  }

  std::vector<float> hit_times[2];  // times of the hits in each [arm]

  // calculate bbc global variables
  if ( _verbose>=10 ) cout << "Hit PMT info " << endl;
  for (int ipmt=0; ipmt<BbcDefs::BBC_N_PMT; ipmt++)
  {
    BbcPmtInfoV1 *bbcpmt = bbcpmts->get_pmt( ipmt );
    int arm = ipmt/64;

    float t_pmt = bbcpmt->get_t();  // hit time of pmt
    float q_pmt = bbcpmt->get_q();  // charge in pmt

    if ( _verbose>=10 ) cout << ipmt << "\t" << t_pmt << endl;
    //if ( m_q[ipmt] > 20 && m_tq[ipmt] > 0. && m_tq[ipmt] < 35. )
    //if ( m_tq[ipmt] > 0. && m_tq[ipmt] < 35. )
    if ( fabs(t_pmt) < 25. )
    {
      hit_times[arm].push_back( t_pmt );
      hevt_bbct[arm]->Fill( t_pmt );

      m_bbcn[arm]++;
      m_bbcq[arm] += q_pmt;

      if ( _verbose>=10 )
      {
        cout << ipmt << "\t" << t_pmt << "\t" << q_pmt << endl;
      }

    }

  }

  if ( _verbose>=10 ) cout << "nhits " << m_bbcn[0] << "\t" << m_bbcn[1] << endl;
  //cout << "bbcte " << m_bbcte[0] << "\t" << m_bbcte[1] << endl;

  for (int iarm = 0; iarm < 2; iarm++)
  {
    if ( hit_times[iarm].empty() )
    {
      //cout << "hit_times size == 0" << endl;
      continue;
    }

    //cout << "EARLIEST " << iarm << endl;
    //cout << "ERROR hit_times size == " << hit_times[iarm].size() << endl;

    std::sort(hit_times[iarm].begin(), hit_times[iarm].end());
    float earliest = hit_times[iarm].at(0);
    //cout << "earliest" << iarm << "\t" << earliest << endl;

    gaussian->SetParameter(0, 5);
    //gaussian->SetParameter(1, earliest);
    //gaussian->SetRange(6, earliest + 5 * 0.05);
    gaussian->SetParameter(1,hevt_bbct[iarm]->GetMean());
    gaussian->SetParameter(2,hevt_bbct[iarm]->GetRMS());
    gaussian->SetRange(hevt_bbct[iarm]->GetMean()-5,hevt_bbct[iarm]->GetMean()+5);

    if ( _verbose ) 
    {
      if ( ac == nullptr )
      {
        ac = new TCanvas("ac","ac",550*1.5,425*1.5);
        ac->Divide(2,1);
      }
      ac->cd(iarm+1);
    }

    hevt_bbct[iarm]->Fit(gaussian, "BNQLR");
    if ( _verbose ) hevt_bbct[iarm]->Draw();

    // m_bbct[iarm] = m_bbct[iarm] / m_bbcn[iarm];
    m_bbct[iarm] = gaussian->GetParameter(1);
    m_bbcte[iarm] = earliest;

    //_bbcout->set_arm(iarm, m_bbcn[iarm], m_bbcq[iarm], m_bbct[iarm]);
  }

  // Get Zvertex, T0
  if (m_bbcn[0] > 0 && m_bbcn[1] > 0)
  {
    // Now calculate zvtx, t0 from best times
    if ( _verbose>=10 ) 
    {
      cout << "Evt " << _eventnum << "\tt0\t" << m_bbct[0] << "\t" << m_bbct[1] << endl;
      cout << "bbcn " << m_bbcn[0] << "\t" << m_bbcn[1] << endl;
      cout << "bbcq " << m_bbcq[0] << "\t" << m_bbcq[1] << endl;
    }
    m_bbcz = (m_bbct[0] - m_bbct[1]) * TMath::C() * 1e-7 / 2.0;   // in cm
    m_bbct0 = (m_bbct[0] + m_bbct[1]) / 2.0;

    // correct z-vertex
    m_bbcz += bz_offset;

    // hard code these for now
    // need study to determine muliplicity dependence
    m_bbczerr = 1.0;    // cm
    m_bbct0err = 0.05;  // ns

    /*
    // Use earliest time
    //cout << "t0\t" << m_bbct[0] << "\t" << m_bbct[1] << endl;
    //cout << "te\t" << m_bbcte[0] << "\t" << m_bbcte[1] << endl;
    m_bbcz = (m_bbcte[0] - m_bbcte[1]) * TMath::C() * 1e-7 / 2.0; // in cm
    m_bbct0 = (m_bbcte[0] + m_bbcte[1]) / 2.0;
    */

    if ( _verbose>10 ) cout << "bbcz " << m_bbcz << endl;
  }

  // Fill rest of BbcOut
  if ( bbcout!=0 )
  {
    for (int iarm=0; iarm<2; iarm++)
    {
      bbcout->set_arm( iarm, get_bbcn(iarm), get_bbcq(iarm), get_bbct(iarm) );
      if (_verbose>10 ) cout <<  get_bbcn(iarm) << "\t" << get_bbcq(iarm) << "\t" << get_bbct(iarm) << endl;
    }

    if ( get_bbcn(0) > 0 && get_bbcn(1) > 0 )
    {
      bbcout->set_t0( get_bbct0(), get_bbct0err() );
      bbcout->set_zvtx( get_bbcz(), get_bbczerr() );
    }

  }

  return 1;
}

// This needs to be reconsidered for 2024 run, hopefully timing instability is fixed by then!
// Only used in online monitoring
int BbcEvent::DoQuickClockOffsetCalib()
{
  for (int ifeech=0; ifeech<256; ifeech++)
  {
    bbcsig[ifeech].SetXY(m_samp[ifeech],m_adc[ifeech]);

    // determine the trig_samp board by board
    int tq = (ifeech/8)%2;   // 0 = T-channel, 1 = Q-channel
    int pmtch = (ifeech/16)*8 + ifeech%8;

    double x, y;
    bbcsig[ifeech].LocMax(x,y);
    h2_tmax[tq]->Fill(x,pmtch);
  }

  if ( h2_tmax[1]->GetEntries() == 128*100 )
  {
    TString name;
    TH1 *h_trigsamp[16]{};
    for (int iboard=0; iboard<16; iboard++)
    {
      name = "h_trigsamp"; name += iboard;
      h_trigsamp[iboard] = h2_tmax[1]->ProjectionX( name, iboard*8 + 1, (iboard+1)*8 );
      int maxbin = h_trigsamp[iboard]->GetMaximumBin();
      TRIG_SAMP[iboard] = h_trigsamp[iboard]->GetBinCenter( maxbin );
      //std::cout << "iboard " << iboard << "\t" << iboard*8+1 << "\t" << (iboard+1)*8 << "\t" << h_trigsamp[iboard]->GetEntries() << std::endl;
      cout << "TRIG_SAMP" << iboard << "\t" << TRIG_SAMP[iboard] << endl;
    }

  }

  return 1;
}


int BbcEvent::Read_Charge_Calib(const char *gainfname)
{
  std::ifstream gainfile( gainfname );

  cout << "Reading gains from " << gainfname << endl;
  int ch;
  float integ, integerr;
  float peak, peakerr;
  float width, widtherr;
  float chi2ndf;
  while ( gainfile >> ch >> integ >> peak >> width >> integerr >> peakerr >> widtherr >> chi2ndf )
  {
    gaincorr[ch] = 1.0/peak;

    //cout << ch << "\t" << peak << endl;
  }

  gainfile.close();

  return 1;
}

// Read in tq t0 offset calibrations
int BbcEvent::Read_TQ_T0_Offsets(const char *t0cal_fname)
{
  ifstream tcalibfile( t0cal_fname );

  cout << "Reading tq_t0 offset calibrations from " << t0cal_fname << endl;

  int pmtnum;
  float meanerr;
  float sigma;
  float sigmaerr;
  for (int ipmt=0; ipmt<BbcDefs::BBC_N_PMT; ipmt++)
  {
    tcalibfile >> pmtnum >> tq_t0_offsets[ipmt] >> meanerr >> sigma >> sigmaerr;
    if ( pmtnum != ipmt )
    {
      cerr << "ERROR, pmtnum != ipmt, " << pmtnum << "\t" << ipmt << endl;
    }
  }

  tcalibfile.close();

  return 1;
}

// Read in tq clk offset calibrations
int BbcEvent::Read_TQ_CLK_Offsets(const char *t0cal_fname)
{
  ifstream tcalibfile( t0cal_fname );

  cout << "Reading tq_clk offset calibrations from " << t0cal_fname << endl;

  int pmtnum;
  for (int ipmt=0; ipmt<BbcDefs::BBC_N_PMT; ipmt++)
  {
    tcalibfile >> pmtnum >> tq_clk_offsets[ipmt];
    if ( pmtnum != ipmt )
    {
      cerr << "ERROR, pmtnum != ipmt, " << pmtnum << "\t" << ipmt << endl;
    }
  }

  tcalibfile.close();

  return 1;
}


// Read in tt clk offset calibrations
int BbcEvent::Read_TT_CLK_Offsets(const char *t0cal_fname)
{
  ifstream tcalibfile( t0cal_fname );

  cout << "Reading tq_clk offset calibrations from " << t0cal_fname << endl;

  int pmtnum;
  for (int ipmt=0; ipmt<BbcDefs::BBC_N_PMT; ipmt++)
  {
    tcalibfile >> pmtnum >> tt_clk_offsets[ipmt];
    if ( pmtnum != ipmt )
    {
      cerr << "ERROR, pmtnum != ipmt, " << pmtnum << "\t" << ipmt << endl;
    }
  }

  tcalibfile.close();

  return 1;
}



