#include "MbdCalibReco.h"
#include "MbdCalib.h"
#include "MbdDefs.h"
#include "MbdPmtContainer.h"
#include "MbdPmtHit.h"
#include "MbdOut.h"

#include <ffamodules/CDBInterface.h>
#include <ffaobjects/EventHeader.h>
#include <ffaobjects/RunHeader.h>
#include <ffarawobjects/Gl1Packet.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TDirectory.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

MbdCalibReco::MbdCalibReco(const std::string &name)
  : SubsysReco(name)
{
}

int MbdCalibReco::Init(PHCompositeNode * /*topNode*/)
{
  _mbdcal = std::make_unique<MbdCalib>();
  _mbdcal->Verbosity( Verbosity() );
  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdCalibReco::InitRun(PHCompositeNode *topNode)
{
  _runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (!_runheader)
  {
    std::cout << PHWHERE << " RunHeader node not found, will use run number 0" << std::endl;
  }

  _runnumber = _runheader ? _runheader->get_RunNumber() : 0;

  getNodes(topNode);

  // Build run directory path and create it
  std::ostringstream oss;
  oss << _caldir << "/" << _runnumber;
  _rundir = oss.str();
  gSystem->Exec(("mkdir -p " + _rundir).c_str());

  if (!_cdbtag.empty())
  {
    // Download baseline calibrations from CDB
    recoConsts::instance()->set_StringFlag("CDB_GLOBALTAG", _cdbtag);
    CDBInterface* cdb = CDBInterface::instance();
    std::string url;

    url = cdb->getUrl("MBD_SAMPMAX");
    if (!url.empty()) { _mbdcal->Download_SampMax(url); }

    url = cdb->getUrl("MBD_PED");
    if (!url.empty()) { _mbdcal->Download_Ped(url); }

    url = cdb->getUrl("MBD_TIMECORR");
    if (!url.empty()) { _mbdcal->Download_TimeCorr(url); }

    url = cdb->getUrl("MBD_SLEWCORR");
    if (!url.empty()) { _mbdcal->Download_SlewCorr(url); }

    std::cout << Name() << ": loaded calibrations from CDB tag " << _cdbtag << std::endl;
  }
  else
  {
    // Load baseline calibrations from local files if they exist
    std::string calfile = _rundir + "/mbd_sampmax.calib";
    if (gSystem->AccessPathName(calfile.c_str()) == 0)
    {
      _mbdcal->Download_SampMax(calfile);
    }
    calfile = _rundir + "/mbd_ped.calib";
    if (gSystem->AccessPathName(calfile.c_str()) == 0)
    {
      _mbdcal->Download_Ped(calfile);
    }

    // Load slew correction for subpass >= 2
    if (_subpass >= 2)
    {
      calfile = _rundir + "/mbd_slewcorr.calib";
      if (gSystem->AccessPathName(calfile.c_str()) == 0)
      {
        _mbdcal->Download_SlewCorr(calfile);
        std::cout << Name() << ": loaded " << calfile << std::endl;
      }
      else
      {
        std::cout << Name() << ": WARNING: " << calfile << " not found" << std::endl;
      }
    }
  }

  // Load t0 offsets for subpass >= 1 (always from local files — outputs of previous subpass)
  if (_subpass >= 1)
  {
    std::string prevpass = "pass" + std::to_string(_subpass - 1) + "_";

    std::string calfile = _rundir + "/" + prevpass + "mbd_tq_t0.calib";
    if (gSystem->AccessPathName(calfile.c_str()) == 0)
    {
      _mbdcal->Download_TQT0(calfile);
      std::cout << Name() << ": loaded " << calfile << std::endl;
    }
    else
    {
      std::cout << Name() << ": WARNING: " << calfile << " not found" << std::endl;
    }

    calfile = _rundir + "/" + prevpass + "mbd_tt_t0.calib";
    if (gSystem->AccessPathName(calfile.c_str()) == 0)
    {
      _mbdcal->Download_TTT0(calfile);
      std::cout << Name() << ": loaded " << calfile << std::endl;
    }
    else
    {
      std::cout << Name() << ": WARNING: " << calfile << " not found" << std::endl;
    }
  }

  // Build bitmask of scaled triggers whose names begin with "MBD N&S"
  _mbias_trigger_mask = 0xfc00;

  // Open output ROOT file
  TDirectory *origdir = gDirectory;

  std::string outfname = _rundir + "/calmbdpass2." + std::to_string(_subpass);
  if (_subpass == 0)
  {
    outfname += "_time-" + std::to_string(_runnumber) + ".root";
  }
  else if (_subpass == 1 || _subpass == 2)
  {
    outfname += "_slew-" + std::to_string(_runnumber) + ".root";
  }
  else
  {
    outfname += "_q-" + std::to_string(_runnumber) + ".root";
  }
  _outfile = std::make_unique<TFile>(outfname.c_str(), "RECREATE");
  if (!_outfile || _outfile->IsZombie())
  {
    std::cerr << PHWHERE << " ERROR: cannot open output file " << outfname << std::endl;
    _outfile.reset();
    return Fun4AllReturnCodes::ABORTRUN;
  }
  std::cout << Name() << ": output file " << outfname << std::endl;

  BookHistograms();

  origdir->cd();

  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdCalibReco::getNodes(PHCompositeNode *topNode)
{
  _evtheader = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if (!_evtheader)
  {
    std::cout << PHWHERE << " EvtHeader not found, will use run number 0" << std::endl;
  }

  _gl1packet = findNode::getClass<Gl1Packet>(topNode,14001);
  if (!_gl1packet)
  {
    _gl1packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    static int counter = 0;
    if ( counter<4 )
    {
      std::cout << PHWHERE << " GL1Packet not found" << std::endl;
    }
  }

  _mbdpmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  if (!_mbdpmts)
  {
    static int counter = 0;
    if ( counter<4 )
    {
      std::cout << PHWHERE << " MbdPmtContainer not found" << std::endl;
    }
  }

  _mbdout = findNode::getClass<MbdOut>(topNode, "MbdOut");
  if (!_mbdout)
  {
    static int counter = 0;
    if ( counter<4 )
    {
      std::cout << PHWHERE << " MbdOut not found" << std::endl;
    }
  }

  _mbdgeom = findNode::getClass<MbdGeom>(topNode, "MbdGeom");
  if (!_mbdgeom)
  {
    static int counter = 0;
    if ( counter<4 )
    {
      std::cout << PHWHERE << " MbdGeom not found" << std::endl;
    }
  }

  if ( !_mbdgeom || !_mbdout || !_mbdpmts || !_gl1packet || !_evtheader )
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void MbdCalibReco::BookHistograms()
{
  // Delete histograms if they have already have been booked.
  if ( h2_tt )
  {
    DeleteHistograms();
    return;
  }

  for (int ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
  {
    std::string sn = std::to_string(ipmt);

    h_tt[ipmt] = new TH1F(("h_tt" + sn).c_str(), ("tt" + sn).c_str(), 7000, -30., 30.);
    h_tt[ipmt]->SetXTitle("ns");

    h_tq[ipmt] = new TH1F(("h_tq" + sn).c_str(), ("tq" + sn).c_str(), 7000, -150., 31. * 17.7623);
    h_tq[ipmt]->SetXTitle("ns");

    h_qp[ipmt] = new TH1F(("h_q" + sn).c_str(), ("q" + sn).c_str(), 3000, -100., 14900.);
    h_qp[ipmt]->SetXTitle("ADC");

    if (_subpass >= 1)
    {
      const int    nbins[2] = {4000, 1100};
      const double xmin[2]  = {-0.5, -5.};
      const double xmax[2]  = {16000. - 0.5, 6.};
      h2_slew[ipmt] = new THnSparseF(("h2_slew" + sn).c_str(), ("slew curve, ch " + sn).c_str(), 2, nbins, xmin, xmax);
      h2_slew[ipmt]->GetAxis(0)->SetTitle("ADC");
      h2_slew[ipmt]->GetAxis(1)->SetTitle("#Delta T (ns)");
    }
    else
    {
      h2_slew[ipmt] = nullptr;
    }
  }

  h2_tt = new TH2F("h2_tt", "ch vs tt", 900, -150., 150., MbdDefs::MBD_N_PMT, -0.5, MbdDefs::MBD_N_PMT - 0.5);
  h2_tt->SetXTitle("tt [ns]");
  h2_tt->SetYTitle("pmt ch");

  h2_tq = new TH2F("h2_tq", "ch vs tq", 900, -150., 150., MbdDefs::MBD_N_PMT, -0.5, MbdDefs::MBD_N_PMT - 0.5);
  h2_tq->SetXTitle("tq [ns]");
  h2_tq->SetYTitle("pmt ch");
}

void MbdCalibReco::DeleteHistograms()
{
  for (int ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
  {
    if ( h_tt[ipmt] )
    {
      delete h_tt[ipmt];
    }
    if ( h_tq[ipmt] )
    {
      delete h_tq[ipmt];
    }
    if ( h_qp[ipmt] )
    {
      delete h_qp[ipmt];
    }
    if ( h2_slew[ipmt] )
    {
      delete h2_slew[ipmt];
    }
  }

  delete h2_tt;
  delete h2_tq;
}

int MbdCalibReco::process_event(PHCompositeNode * /*topNode*/)
{
  // Require a scaled "MBD N&S" trigger
  if (_mbias_trigger_mask != 0)
  {
    uint64_t strig = _gl1packet->getScaledVector();  // scaled trigger only
    if ( (strig&_mbias_trigger_mask)==0 )
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  std::array<Float_t, MbdDefs::MBD_N_ARMS> armtime{};
  armtime.fill(0);
  std::array<Float_t, MbdDefs::MBD_N_ARMS> nhit{};
  nhit.fill(0);

  Float_t zvtx = _mbdout->get_zvtx();
  // Vertex cut for subpass >= 1
  if ( _subpass >= 1 )
  {
    if (std::abs(zvtx) > 60.)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  for (int iarm=0; iarm<2; iarm++)
  {
    armtime[iarm] = _mbdout->get_time(iarm);
    nhit[iarm] = _mbdout->get_npmt(iarm);
  }

  for (int ipmt=0; ipmt < _mbdpmts->get_npmt(); ipmt++)
  {
    MbdPmtHit *pmt = _mbdpmts->get_pmt(ipmt);
    if ( !pmt )
    {
      continue;
    }

    Short_t pmtno = pmt->get_pmt();
    if ( pmtno<0 || pmtno>=MbdDefs::MBD_N_PMT )
    {
      static int counter = 0;
      if ( counter<10 )
      {
        std::cerr << PHWHERE << " invalide pmt no " << pmtno << std::endl;
        counter++;
      }
      continue;
    }

    Float_t q = pmt->get_q();
    Float_t tt  = pmt->get_tt();
    Float_t tq  = pmt->get_tq();

    h_tt[pmtno]->Fill( tt );
    h2_tt->Fill( tt, pmtno );
    h_tq[pmtno]->Fill( tq );
    h2_tq->Fill( tq, pmtno );

    // Fill charge histogram for in-time hits
    if ( std::abs(tt)<26.0 && q > 0.)
    {
      h_qp[pmtno]->Fill( q );
    }

    int arm = _mbdgeom->get_arm( pmtno );

    // Fill slew histogram for subpass >= 1
    if (_subpass >= 1 && h2_slew[pmtno])
    {
      if (nhit[arm] >= 2. && q > 0.)
      {
        float dt = tt - armtime[arm];
        const double coords[2] = {q, dt};
        h2_slew[pmtno]->Fill(coords);
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdCalibReco::EndRun(const int /*runnumber*/)
{
  if (!_outfile)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // Write histograms to output file
  _outfile->cd();
  h2_tt->Write();
  h2_tq->Write();
  for (int ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
  {
    h_tt[ipmt]->Write();
    h_tq[ipmt]->Write();
    h_qp[ipmt]->Write();
    if (h2_slew[ipmt])
    {
      h2_slew[ipmt]->Write();
    }
  }

  _outfile->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

