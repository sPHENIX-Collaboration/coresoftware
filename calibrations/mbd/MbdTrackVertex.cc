#include "MbdTrackVertex.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <ffarawobjects/Gl1Packet.h>

#include <ffaobjects/EventHeader.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

constexpr GlobalVertex::VTXTYPE trkType = GlobalVertex::SVTX;
constexpr GlobalVertex::VTXTYPE mbdType = GlobalVertex::MBD;
//____________________________________________________________________________..
MbdTrackVertex::MbdTrackVertex(const std::string &name):
 SubsysReco(name)
{
  std::cout << "MbdTrackVertex::MbdTrackVertex(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
MbdTrackVertex::~MbdTrackVertex()
{
  std::cout << "MbdTrackVertex::~MbdTrackVertex() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int MbdTrackVertex::Init(PHCompositeNode * /*topNode*/)
{
  outFile = new TFile(outFileName.c_str(), "RECREATE");
  if (_treeflag)
  {
    outTree = new TTree("mz", "MBD-TRK ZVTX");
    outTree->OptimizeBaskets();
    outTree->SetAutoSave(-5e6);

    outTree->Branch("evt", &_evt, "evt/I");
    outTree->Branch("mbdz", &_mbdVertex, "mbdz/F");
    outTree->Branch("trkz", &_trackerVertex, "trkz/F");
    outTree->Branch("ntrks", &_nTracks, "ntrks/i");
    outTree->Branch("nbz", &_nMBDVertex, "nbz/i");
    outTree->Branch("ntz", &_nTRKVertex, "ntz/i");
  }

  // h_mbdtrkz: dz = _mbdVertex - trackerVertex, range (-15,15) cm, 0.25 cm bins → 120 bins
  h_mbdtrkz = new TH1F("h_mbdtrkz", "MBD - Tracker z-vertex;dz (cm);Counts", 120, -15., 15.);
  h_bz   = new TH1F("h_bz",   "MBD z-vertex;z (cm);Counts", 400, -20., 20.);
  h_trkz = new TH1F("h_trkz", "Tracker z-vertex;z (cm);Counts", 400, -20., 20.);

  // h2_mbdtrkz: THnSparseF, x = _trackerVertex, y = _mbdVertex, (-20,20) cm, 0.1 cm bins → 400 bins each
  const int    nbins2[2] = {400, 400};
  const double xmin2[2]  = {-20., -20.};
  const double xmax2[2]  = { 20.,  20.};
  h2_mbdtrkz = new THnSparseF("h2_mbdtrkz", "MBD vs Tracker z-vertex", 2, nbins2, xmin2, xmax2);
  h2_mbdtrkz->GetAxis(0)->SetTitle("Tracker z (cm)");
  h2_mbdtrkz->GetAxis(1)->SetTitle("MBD z (cm)");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MbdTrackVertex::process_event(PHCompositeNode *topNode)
{
  EventHeader *eventheader = findNode::getClass<EventHeader>(topNode, "EventHeader");
  _evt = eventheader ? eventheader->get_EvtSequence() : -1;

  if (_gl1_trigmask != 0)
  {
    Gl1Packet *gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (gl1)
    {
      if ((gl1->getScaledVector() & _gl1_trigmask) == 0)
      {
        return Fun4AllReturnCodes::DISCARDEVENT;
      }
    }
    else
    {
      std::cout << PHWHERE << " GL1Packet node not found; discarding event because trigger masking was requested" << std::endl;
      return Fun4AllReturnCodes::DISCARDEVENT;
    }
  }

  MbdVertexMap *m_dst_mbdvertexmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
  SvtxVertexMap *m_dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");

  GlobalVertexMap *globalvertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!m_dst_mbdvertexmap || !m_dst_vertexmap || !globalvertexmap)
  {
    std::cout << PHWHERE << " missing required vertex node(s)" << std::endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  _mbdVertex = _trackerVertex = std::numeric_limits<float>::quiet_NaN();
  _nTracks = _nMBDVertex = _nTRKVertex = std::numeric_limits<unsigned int>::quiet_NaN();

  _hasMBD = false;
  _hasTRK = false;

  for (GlobalVertexMap::ConstIter iter = globalvertexmap->begin(); iter != globalvertexmap->end(); ++iter)
  {
    GlobalVertex *gvertex = iter->second;

    if (gvertex->count_vtxs(mbdType) != 0)
    {
      _hasMBD = true;

      auto mbditer = gvertex->find_vertexes(mbdType);
      auto mbdvertexvector = mbditer->second;

      _nMBDVertex = mbdvertexvector.size();
      for (auto &vertex : mbdvertexvector)
      {
        MbdVertex *m_dst_vertex = m_dst_mbdvertexmap->find(vertex->get_id())->second;
        _mbdVertex = m_dst_vertex->get_z();
      }
    }

    if (gvertex->count_vtxs(trkType) != 0)
    { 
      _hasTRK = true; 

      auto trkiter = gvertex->find_vertexes(trkType);
      auto trkvertexvector = trkiter->second;

      _nTRKVertex = trkvertexvector.size();
      for (auto &vertex : trkvertexvector)
      {
        SvtxVertex *m_dst_vertex = m_dst_vertexmap->find(vertex->get_id())->second;
        if ( m_dst_vertex->get_beam_crossing() != 0 ) 
        {
          continue;
        }
        if ( m_dst_vertex->size_tracks() > _nTracks)
        {
          _trackerVertex = m_dst_vertex->get_z();
          _nTracks = m_dst_vertex->size_tracks();
        }
        if (_nTracks == 0)
        {
          _hasTRK = false; 
        }
      }
    }
  }

  if (_hasMBD)
  {
    h_bz->Fill(_mbdVertex);
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
  outFile->Write();
  outFile->Close();
  delete outFile;

  return Fun4AllReturnCodes::EVENT_OK;
}

