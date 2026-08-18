#include "PHG4InEventToTruthInfo.h"

#include "PHG4InEvent.h"
#include "PHG4MCProcessDefs.h"
#include "PHG4Particle.h"
#include "PHG4Particlev2.h"
#include "PHG4Particlev3.h"
#include "PHG4Showerv1.h"
#include "PHG4TruthInfoContainer.h"
#include "PHG4VtxPoint.h"
#include "PHG4VtxPointv2.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TDatabasePDG.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>

PHG4InEventToTruthInfo::PHG4InEventToTruthInfo(const std::string &name)
  : SubsysReco(name)
{
}

int PHG4InEventToTruthInfo::InitRun(PHCompositeNode *topNode)
{
  CreateNodeTree(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4InEventToTruthInfo::process_event(PHCompositeNode *topNode)
{
  PHG4InEvent *inevent = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!inevent)
  {
    if (Verbosity() > 0)
    {
      std::cout << PHWHERE << "PHG4INEVENT node not found" << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  std::map<int, int> input_to_truth_vtxid;
  const auto vtxrange = inevent->GetVertices();
  for (auto vtxiter = vtxrange.first; vtxiter != vtxrange.second; ++vtxiter)
  {
    const int input_vtxid = vtxiter->first;
    const int truth_vtxid = truthinfo->maxvtxindex() + 1;
    PHG4VtxPoint *truth_vtx = dynamic_cast<PHG4VtxPoint *>(vtxiter->second->CloneMe());
    truth_vtx->set_id(truth_vtxid);

    const auto inserted = truthinfo->AddVertex(truth_vtxid, truth_vtx);
    if (inserted == truthinfo->GetVtxRange().second)
    {
      std::cout << PHWHERE << " Failure to insert vertex " << truth_vtxid << " into TruthInfo" << std::endl;
      truth_vtx->identify();
      gSystem->Exit(1);
      exit(1);
    }

    input_to_truth_vtxid[input_vtxid] = truth_vtxid;
  }

  const auto particlerange = inevent->GetParticles();
  // reverse to match how G4 inserts particles into the truth table (last particle first)
  for (auto particleiter = std::make_reverse_iterator(particlerange.second);
       particleiter != std::make_reverse_iterator(particlerange.first);
       ++particleiter)
  {
    const int input_vtxid = particleiter->first;
    PHG4Particle *input_particle = particleiter->second;
    if (!input_particle)
    {
      continue;
    }

    auto vtxid_iter = input_to_truth_vtxid.find(input_vtxid);
    if (vtxid_iter == input_to_truth_vtxid.end())
    {
      std::cout << PHWHERE << "no vertex " << input_vtxid
                << " found for PHG4InEvent particle" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    const int truth_trackid = truthinfo->maxtrkindex() + 1;
    const int truth_vtxid = vtxid_iter->second;

    PHG4Particle *truth_particle = makeParticleCopy(input_particle);
    truth_particle->set_track_id(truth_trackid);
    truth_particle->set_vtx_id(truth_vtxid);
    truth_particle->set_parent_id(0);
    truth_particle->set_primary_id(truth_trackid);

    const auto inserted = truthinfo->AddParticle(truth_trackid, truth_particle);
    if (inserted == truthinfo->GetParticleRange().second)
    {
      std::cout << PHWHERE << " Failure to insert particle " << truth_trackid << " into TruthInfo" << std::endl;
      truth_particle->identify();
      gSystem->Exit(1);
      exit(1);
    }

    const int embed_flag = inevent->isEmbeded(input_particle);
    if (embed_flag)
    {
      truthinfo->AddEmbededTrkId(truth_trackid, embed_flag);
      truthinfo->AddEmbededVtxId(truth_vtxid, embed_flag);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

PHG4Particle *PHG4InEventToTruthInfo::makeParticleCopy(PHG4Particle *particle)
{
  if (particle->get_name().empty())
  {
    particle->set_name(TDatabasePDG::Instance()->GetParticle(particle->get_pid())->GetName());
  }
  PHG4Particle *particle_return{nullptr};
  if (particle->isIon())
  {
    particle_return = new PHG4Particlev3(particle);
  }
  else
  {
    particle_return = new PHG4Particlev2(particle);
  }

  if (!std::isfinite(particle_return->get_e()))
  {
    double m = TDatabasePDG::Instance()->GetParticle(particle->get_pid())->Mass();
    double e = sqrt(particle->get_px() * particle->get_px() + particle->get_py() * particle->get_py() + particle->get_pz() * particle->get_pz() + m * m);
    particle_return->set_e(e);
  }
  return particle_return;
}

int PHG4InEventToTruthInfo::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST node not found" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(dstNode, "G4TruthInfo");
  if (!truthinfo)
  {
    truthinfo = new PHG4TruthInfoContainer();
    dstNode->addNode(new PHIODataNode<PHObject>(truthinfo, "G4TruthInfo", "PHObject"));
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
