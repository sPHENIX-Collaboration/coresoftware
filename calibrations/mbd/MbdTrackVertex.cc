#include "MbdTrackVertex.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <ffarawobjects/Gl1Packet.h>

#include <ffaobjects/EventHeader.h>
#include <ffaobjects/RunHeader.h>

#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <mbd/MbdOut.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TTree.h>

#include <cmath>
#include <iostream>

//____________________________________________________________________________..
MbdTrackVertex::MbdTrackVertex(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int MbdTrackVertex::Init(PHCompositeNode * /*topNode*/)
{
  if ( outFile != nullptr )
  {
    outFile->Close();
    delete outFile;
  }

  outFile = new TFile(outFileName.c_str(), "RECREATE");
  if (_treeflag)
  {
    outTree = new TTree("mz", "MBD-TRK ZVTX");
    outTree->OptimizeBaskets();
    outTree->SetAutoSave(-5e6);

    outTree->Branch("run", &_run, "run/I");
    outTree->Branch("evt", &_evt, "evt/I");
    outTree->Branch("strig", &_strig, "strig/l");
    outTree->Branch("mbdz", &_mbdVertex, "mbdz/F");
    outTree->Branch("mbdt0", &_mbdT0, "mbdt0/F");
    outTree->Branch("trkz", &_trackerVertex, "trkz/F");
    outTree->Branch("trkx", &_trackerX, "trkx/F");
    outTree->Branch("trky", &_trackerY, "trky/F");
    outTree->Branch("trkchi2ndf", &_trackerChi2ndf, "trkchi2ndf/F");
    outTree->Branch("nbz", &_nMBDVertex, "nbz/i");
    outTree->Branch("ntz", &_nTRKVertex, "ntz/i");
    outTree->Branch("ntzall", &_nTRKVertexAll, "ntzall/i");
    outTree->Branch("ntrks", &_nTracks, "ntrks/i");
    outTree->Branch("bns", &_nMBD[0], "bns/S");
    outTree->Branch("bnn", &_nMBD[1], "bnn/S");
    outTree->Branch("bqs", &_qMBD[0], "bqs/F");
    outTree->Branch("bqn", &_qMBD[1], "bqn/F");
  }

  // h_mbdtrkz: dz = _mbdVertex - trackerVertex, range (-15,15) cm, 0.25 cm bins → 120 bins
  h_mbdtrkz = new TH1F("h_mbdtrkz", "MBD - Tracker z-vertex;dz (cm);Counts", 120, -15., 15.);
  h_bz = new TH1F("h_bz", "MBD z-vertex;z (cm);Counts", 400, -20., 20.);
  h_trkz = new TH1F("h_trkz", "Tracker z-vertex;z (cm);Counts", 400, -20., 20.);
  // h_bt0: MBD t0, range (-10,10) ns, 0.05 ns bins → 400 bins
  h_bt0 = new TH1F("h_bt0", "MBD t0;t0 (ns);Counts", 400, -10., 10.);

  // h2_mbdtrkz: THnSparseF, x = _trackerVertex, y = _mbdVertex, (-20,20) cm, 0.1 cm bins → 400 bins each
  const int nbins2[2] = {400, 400};
  const double xmin2[2] = {-20., -20.};
  const double xmax2[2] = {20., 20.};
  h2_mbdtrkz = new THnSparseF("h2_mbdtrkz", "MBD vs Tracker z-vertex", 2, nbins2, xmin2, xmax2);
  h2_mbdtrkz->GetAxis(0)->SetTitle("Tracker z (cm)");
  h2_mbdtrkz->GetAxis(1)->SetTitle("MBD z (cm)");

  // h2_dzvstrkz: x = _trackerVertex (-20,20) cm, 0.5 cm bins; y = dz (-25,25) cm, 0.25 cm bins
  h2_dzvstrkz = new TH2F("h2_dzvstrkz", "MBD - Tracker z-vertex vs Tracker z;Tracker z (cm);MBD z - Tracker z (cm)",
                         80, -20., 20., 200, -25., 25.);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MbdTrackVertex::process_event(PHCompositeNode *topNode)
{
  RunHeader *runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  _run = runheader ? runheader->get_RunNumber() : -1;

  EventHeader *eventheader = findNode::getClass<EventHeader>(topNode, "EventHeader");
  _evt = eventheader ? eventheader->get_EvtSequence() : -1;

  if (Verbosity())
  {
    std::cout << "evt " << _evt << std::endl;
  }

  // GL1 is packet 14001 in DST_CALOFITTING, GL1RAWHIT in DST_TRKR_TRACKS
  Gl1Packet *gl1 = findNode::getClass<Gl1Packet>(topNode, 14001);
  if (!gl1)
  {
    gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
  }
  if (!gl1)
  {
    gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
  }
  _strig = gl1 ? gl1->getScaledVector() : 0;

  if (_gl1_trigmask != 0)
  {
    if (!gl1)
    {
      std::cout << PHWHERE << " Gl1Packet not found; discarding event because trigger masking was requested" << std::endl;
      return Fun4AllReturnCodes::DISCARDEVENT;
    }
    if ((_strig & _gl1_trigmask) == 0)
    {
      return Fun4AllReturnCodes::DISCARDEVENT;
    }
  }

  MbdOut *mbdout = findNode::getClass<MbdOut>(topNode, "MbdOut");
  MbdVertexMap *mbdvertexmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
  SvtxVertexMap *svtxvertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");

  if (!mbdvertexmap || !svtxvertexmap)
  {
    std::cout << PHWHERE << " missing MbdVertexMap or SvtxVertexMap" << std::endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  _mbdVertex = _mbdT0 = std::numeric_limits<float>::quiet_NaN();
  _trackerVertex = _trackerX = _trackerY = _trackerChi2ndf = std::numeric_limits<float>::quiet_NaN();
  _nTracks = 0;
  _nMBDVertex = 0;
  _nTRKVertex = 0;
  _nTRKVertexAll = 0;

  if (mbdout)
  {
    _nMBD[0] = mbdout->get_npmt(0);
    _nMBD[1] = mbdout->get_npmt(1);
    _qMBD[0] = mbdout->get_q(0);
    _qMBD[1] = mbdout->get_q(1);
  }
  else
  {
    _nMBD[0] = 0;
    _nMBD[1] = 0;
    _qMBD[0] = std::numeric_limits<Float_t>::quiet_NaN();
    _qMBD[1] = std::numeric_limits<Float_t>::quiet_NaN();
  }

  // MBD vertex (MbdReco only fills a valid z when both arms fire)
  _hasMBD = false;
  _nMBDVertex = mbdvertexmap->size();
  for (auto &[key, mbdvtx] : *mbdvertexmap)
  {
    if (std::isfinite(mbdvtx->get_z()))
    {
      _mbdVertex = mbdvtx->get_z();
      _mbdT0 = mbdvtx->get_t();
      _hasMBD = true;
      break;
    }
  }

  // Tracker vertex: the one with the most tracks in the selected crossing
  _hasTRK = false;
  _nTRKVertexAll = svtxvertexmap->size();
  for (auto &[key, svtxvtx] : *svtxvertexmap)
  {
    if (svtxvtx->get_beam_crossing() != _trk_crossing)
    {
      continue;
    }
    ++_nTRKVertex;
    const unsigned int ntrk = svtxvtx->size_tracks();
    if (ntrk > _nTracks)
    {
      _nTracks = ntrk;
      _trackerVertex = svtxvtx->get_z();
      _trackerX = svtxvtx->get_x();
      _trackerY = svtxvtx->get_y();
      _trackerChi2ndf = (svtxvtx->get_ndof() > 0) ? svtxvtx->get_chisq() / svtxvtx->get_ndof() : std::numeric_limits<float>::quiet_NaN();
    }
  }
  if (_nTracks >= _min_ntracks && std::isfinite(_trackerVertex))
  {
    _hasTRK = true;
  }

  if (Verbosity())
  {
    std::cout << "  mbdz " << _mbdVertex << " trkz " << _trackerVertex
              << " ntrks " << _nTracks << " ntz " << _nTRKVertex << "/" << _nTRKVertexAll << std::endl;
  }

  if (_hasMBD)
  {
    h_bz->Fill(_mbdVertex);
    if (std::isfinite(_mbdT0))
    {
      h_bt0->Fill(_mbdT0);
    }
  }
  if (_hasTRK)
  {
    h_trkz->Fill(_trackerVertex);
  }

  if (_hasMBD && _hasTRK)
  {
    h_mbdtrkz->Fill(_mbdVertex - _trackerVertex);
    const double coords[2] = {_trackerVertex, _mbdVertex};
    h2_mbdtrkz->Fill(coords);
    h2_dzvstrkz->Fill(_trackerVertex, _mbdVertex - _trackerVertex);
  }

  if (_treeflag)
  {
    outTree->Fill();
  }

  ++_counter;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MbdTrackVertex::End(PHCompositeNode * /*topNode*/)
{
  std::cout << Name() << " processed " << _counter << " events" << std::endl;
  outFile->Write();
  outFile->Close();
  delete outFile;

  return Fun4AllReturnCodes::EVENT_OK;
}
