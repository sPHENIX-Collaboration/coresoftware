#include "G4VtxNtuple.h"

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>

#include <TNtuple.h>

#include <sstream>

using namespace std;

G4VtxNtuple::G4VtxNtuple(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , m_FileName(filename)
{
}

G4VtxNtuple::~G4VtxNtuple()
{
  delete hm;
}

int G4VtxNtuple::Init(PHCompositeNode *)
{
  hm = new Fun4AllHistoManager(Name());
  ntup = new TNtuple("vtxntup", "G4Vtxs", "vx:vy:vz");
  hm->registerHisto(ntup);
  return 0;
}

int G4VtxNtuple::process_event(PHCompositeNode *topNode)
{
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (truthinfo)
  {
    PHG4VtxPoint *gvertex = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
    ntup->Fill(gvertex->get_x(), gvertex->get_y(), gvertex->get_z());
  }
  return 0;
}

int G4VtxNtuple::End(PHCompositeNode */*topNode*/)
{
  hm->dumpHistos(m_FileName, "RECREATE");
  return 0;
}
