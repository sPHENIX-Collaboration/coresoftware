#include "TruthNeutralMesonBuilder.h"

#include "TruthNeutralMesonContainer.h"
#include "TruthNeutralMesonv1.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#pragma GCC diagnostic pop
#include <HepMC/GenParticle.h>

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

namespace
{
  inline float safe_pt(float px, float py)
  {
    return std::sqrt(px * px + py * py);
  }
  inline float safe_p(float px, float py, float pz)
  {
    return std::sqrt(px * px + py * py + pz * pz);
  }
  inline float safe_phi(float px, float py)
  {
    return std::atan2(py, px);
  }
  inline float safe_eta(float p, float pz)
  {
    const float denom = p - std::fabs(pz);
    if (p <= 0 || denom <= 0)
    {
      return 0;
    }
    const float num = p + std::fabs(pz);
    return 0.5 * std::log(num / denom) * (pz < 0 ? -1 : 1);
  }
}  // namespace

TruthNeutralMesonBuilder::TruthNeutralMesonBuilder(const std::string &name)
  : SubsysReco(name)
{
}

int TruthNeutralMesonBuilder::Init(PHCompositeNode *topNode)
{
  CreateNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int TruthNeutralMesonBuilder::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << " DST node missing." << std::endl;
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
    std::cout << "Added PHComposite node: DST" << std::endl;
  }

  _container = findNode::getClass<TruthNeutralMesonContainer>(dstNode, m_truthnmeson_node_name);
  if (!_container)
  {
    _container = new TruthNeutralMesonContainer();
    PHIODataNode<PHObject> *node = new PHIODataNode<PHObject>(_container, m_truthnmeson_node_name, "PHObject");
    dstNode->addNode(node);
    std::cout << "Added TruthNeutralMesonContainer node" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TruthNeutralMesonBuilder::process_event(PHCompositeNode *topNode)
{
  genevtmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!genevtmap && !truthinfo)
  {
    std::cout << PHWHERE << " Missing both PHHepMCGenEventMap or G4TruthInfo. Abort run" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _container = findNode::getClass<TruthNeutralMesonContainer>(topNode, m_truthnmeson_node_name);
  if (!_container)
  {
    std::cout << PHWHERE << " Missing output node " << m_truthnmeson_node_name << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_seen_mother_keys.clear();
  hepmc_by_barcode.clear();
  m_seen_barcodes.clear();

  if (genevtmap)
  {
    for (auto &it : *genevtmap)
    {
      PHHepMCGenEvent *g = it.second;
      if (!g || !g->getEvent())
      {
        continue;
      }
      const int embedding_id = g->get_embedding_id();
      HepMC::GenEvent *evt = g->getEvent();
      for (HepMC::GenEvent::particle_const_iterator hepmcgenit = evt->particles_begin(); hepmcgenit != evt->particles_end(); ++hepmcgenit)
      {
        HepMC::GenParticle *p = *hepmcgenit;
        hepmc_by_barcode[p->barcode()] = std::make_pair(p, embedding_id);
      }
    }
  }

  if (truthinfo)
  {
    g4_by_barcode.clear();
    g4_by_id.clear();
    g4_children.clear();

    PHG4TruthInfoContainer::ConstRange range = truthinfo->GetParticleRange();
    for (auto it = range.first; it != range.second; ++it)
    {
      PHG4Particle *p = it->second;
      if (!p)
      {
        continue;
      }
      g4_by_id[p->get_track_id()] = p;
      g4_by_barcode[p->get_barcode()] = p;
      g4_children[p->get_parent_id()].push_back(p);
    }
  }

  auto fill_kin_from_hepmc = [](HepMC::GenParticle *p, float &e, float &pt, float &eta, float &phi, float &magp)
  {
    const auto &m = p->momentum();
    const float px = m.px();
    const float py = m.py();
    const float pz = m.pz();
    const float ee = m.e();
    const float pp = safe_p(px, py, pz);
    e = ee;
    pt = safe_pt(px, py);
    phi = safe_phi(px, py);
    eta = safe_eta(pp, pz);
    magp = pp;
  };

  auto fill_kin_from_g4 = [](PHG4Particle *p, float &e, float &pt, float &eta, float &phi, float &magp)
  {
    const float px = p->get_px();
    const float py = p->get_py();
    const float pz = p->get_pz();
    const float ee = p->get_e();
    const float pp = safe_p(px, py, pz);
    e = ee;
    pt = safe_pt(px, py);
    phi = safe_phi(px, py);
    eta = safe_eta(pp, pz);
    magp = pp;
  };

  auto is_hadron = [](int id)
  {
    int a = std::abs(id);
    return (a > 100 && a < 100000);
  };

  auto classify_origin = [&](int pid, int barcode, bool &is_prompt, int &parent_pid, bool &pi0_has_eta)
  {
    is_prompt = true;
    parent_pid = 0;
    pi0_has_eta = false;

    HepMC::GenParticle *cur = nullptr;
    auto itH = hepmc_by_barcode.find(barcode);
    if (itH != hepmc_by_barcode.end())
    {
      cur = itH->second.first;
    }

    while (cur)
    {
      HepMC::GenVertex *pv = cur->production_vertex();
      if (!pv)
      {
        break;
      }

      if (pv->particles_in_const_begin() == pv->particles_in_const_end())
      {
        break;
      }

      HepMC::GenParticle *par = *(pv->particles_in_const_begin());
      if (!par)
      {
        break;
      }

      int tmp_parentid = par->pdg_id();
      if (!is_hadron(tmp_parentid))
      {
        break;
      }

      if (is_prompt)
      {
        parent_pid = par->pdg_id();
      }
      is_prompt = false;
      if (pid == 111 && std::abs(parent_pid) == 221)
      {
        pi0_has_eta = true;
      }
      cur = par;
    }

    if (!is_prompt)
    {
      return;
    }

    PHG4Particle *g4 = nullptr;
    auto itG = g4_by_barcode.find(barcode);
    if (itG != g4_by_barcode.end())
    {
      g4 = itG->second;
    }

    while (g4 && g4->get_parent_id() > 0)
    {
      PHG4Particle *par = nullptr;
      auto itP = g4_by_id.find(g4->get_parent_id());
      if (itP == g4_by_id.end())
      {
        break;
      }
      par = itP->second;
      if (is_prompt)
      {
        parent_pid = par->get_pid();
      }
      is_prompt = false;
      if (pid == 111 && std::abs(parent_pid) == 221)
      {
        pi0_has_eta = true;
      }
      g4 = par;
    }
  };

  auto collect_photons_hepmc = [&](HepMC::GenParticle *mom, std::vector<HepMC::GenParticle *> &gammas)
  {
    gammas.clear();
    HepMC::GenVertex *vtx = mom->end_vertex();
    if (!vtx)
    {
      return;
    }
    for (HepMC::GenVertex::particles_out_const_iterator it = vtx->particles_out_const_begin(); it != vtx->particles_out_const_end(); ++it)
    {
      HepMC::GenParticle *ch = *it;
      if (!ch)
      {
        continue;
      }
      if (std::abs(ch->pdg_id()) == 22)
      {
        gammas.push_back(ch);
      }
    }
  };

  auto collect_photons_g4 = [&](int mother_barcode, std::vector<PHG4Particle *> &gammas)
  {
    gammas.clear();
    auto itMom = g4_by_barcode.find(mother_barcode);
    if (itMom == g4_by_barcode.end())
    {
      return;
    }
    PHG4Particle *g4mom = itMom->second;
    auto itKids = g4_children.find(g4mom->get_track_id());
    if (itKids == g4_children.end())
    {
      return;
    }
    for (auto *c : itKids->second)
    {
      if (!c)
      {
        continue;
      }
      if (std::abs(c->get_pid()) == 22)
      {
        gammas.push_back(c);
      }
    }
  };

  if (genevtmap)
  {
    for (auto &it : *genevtmap)
    {
      PHHepMCGenEvent *g = it.second;
      if (!g || !g->getEvent())
      {
        continue;
      }
      const int embedding = g->get_embedding_id();
      HepMC::GenEvent *evt = g->getEvent();

      for (HepMC::GenEvent::particle_const_iterator hepmcgenit = evt->particles_begin(); hepmcgenit != evt->particles_end(); ++hepmcgenit)
      {
        HepMC::GenParticle *p = *hepmcgenit;
        if (!p)
        {
          continue;
        }
        const int pid = p->pdg_id();
        if ((pid != 111 && pid != 221) || (pid == 111 && !m_save_pi0) || (pid == 221 && !m_save_eta))
        {
          continue;
        }

        const int bc = p->barcode();
        const std::pair<int, int> key(embedding, bc);
        if (m_seen_mother_keys.contains(key))
        {
          continue;
        }
        m_seen_mother_keys.insert(key);
        m_seen_barcodes.insert(bc);

        auto *mes = new TruthNeutralMesonv1();
        mes->set_pid(pid);

        float e;
        float pt;
        float eta;
        float phi;
        float magp;
        fill_kin_from_hepmc(p, e, pt, eta, phi, magp);
        mes->set_e(e);
        mes->set_pt(pt);
        mes->set_eta(eta);
        mes->set_phi(phi);
        mes->set_p(magp);

        bool is_prompt = false;
        bool pi0_has_eta = false;
        int parent_pid = 0;
        classify_origin(pid, bc, is_prompt, parent_pid, pi0_has_eta);
        mes->set_prompt(is_prompt);
        mes->set_parent_pid(parent_pid);
        if (pid == 111)
        {
          mes->set_mother_is_eta(pi0_has_eta);
        }
        else if (pid == 221)
        {
          bool eta_3pi0 = false;
          bool eta_pi0pipm = false;
          classify_eta_decay(bc, p, eta_3pi0, eta_pi0pipm);
          mes->set_eta_3pi0(eta_3pi0);
          mes->set_eta_pi0pipm(eta_pi0pipm);
        }

        std::vector<HepMC::GenParticle *> gamsH;
        collect_photons_hepmc(p, gamsH);

        if (gamsH.size() == 2)
        {
          for (auto *gph : gamsH)
          {
            float ge;
            float gpt;
            float geta;
            float gphi;
            float gp;
            fill_kin_from_hepmc(gph, ge, gpt, geta, gphi, gp);

            bool isconverted = false;
            auto itg = g4_by_barcode.find(gph->barcode());
            if (itg != g4_by_barcode.end() && itg->second)
            {
              int gphtrackid = itg->second->get_track_id();
              isconverted = FindConversion(truthinfo, gphtrackid, ge);
            }
            mes->add_photon(ge, gpt, geta, gphi, gp, isconverted);
          }
        }
        else if (truthinfo)
        {
          std::vector<PHG4Particle *> gamsG4;
          collect_photons_g4(bc, gamsG4);
          if (Verbosity() > 0)
          {
            if (gamsG4.size() > 2)
            {
              std::cout << " In main loop / Geant part more than 2 photon check..." << std::endl;
              for (auto *gph : gamsG4)
              {
                std::cout << " mother pid " << pid << " daughter pid " << gph->get_pid() << ", parent id " << gph->get_parent_id() << std::endl;
              }
            }
          }

          if (gamsG4.size() == 2)
          {
            for (auto *gph : gamsG4)
            {
              float ge;
              float gpt;
              float geta;
              float gphi;
              float gp;
              fill_kin_from_g4(gph, ge, gpt, geta, gphi, gp);
              bool isconverted = false;
              int gphtrackid = gph->get_track_id();
              isconverted = FindConversion(truthinfo, gphtrackid, ge);
              mes->add_photon(ge, gpt, geta, gphi, gp, isconverted);
            }
          }
        }
        _container->AddMeson(mes);
      }
    }
  }

  if (truthinfo)
  {
    PHG4TruthInfoContainer::ConstRange range = truthinfo->GetParticleRange();
    for (auto it = range.first; it != range.second; ++it)
    {
      PHG4Particle *p = it->second;
      if (!p)
      {
        continue;
      }
      const int pid = p->get_pid();
      if ((pid != 111 && pid != 221) || (pid == 111 && !m_save_pi0) || (pid == 221 && !m_save_eta))
      {
        continue;
      }

      const int bc = p->get_barcode();
      if (m_seen_barcodes.contains(bc))
      {
        continue;
      }
      m_seen_barcodes.insert(bc);
      const std::pair<int, int> key(0, bc);
      if (m_seen_mother_keys.contains(key))
      {
        continue;
      }
      m_seen_mother_keys.insert(key);

      auto *mes = new TruthNeutralMesonv1();
      mes->set_pid(pid);

      float e;
      float pt;
      float eta;
      float phi;
      float magp;
      fill_kin_from_g4(p, e, pt, eta, phi, magp);
      mes->set_e(e);
      mes->set_pt(pt);
      mes->set_eta(eta);
      mes->set_phi(phi);
      mes->set_p(magp);

      bool is_prompt = false;
      bool pi0_has_eta = false;
      int parent_pid = 0;
      classify_origin(pid, bc, is_prompt, parent_pid, pi0_has_eta);
      if (parent_pid == pid || parent_pid == 22)
      {
        continue;
      }

      auto g4parents = g4_by_id.find(p->get_parent_id());
      if (g4parents == g4_by_id.end())
      {
        continue;
      }
      
      auto *gp = g4parents->second;
      if (!gp)
      {
        continue;
      }
      if (gp->get_track_id() < 0)
      {
        continue;
      }
      int parentbc = gp->get_barcode();

      auto itH = hepmc_by_barcode.find(parentbc);
      if (itH == hepmc_by_barcode.end())
      {
        continue;
      }

      mes->set_prompt(is_prompt);
      mes->set_parent_pid(parent_pid);
      if (pid == 111)
      {
        mes->set_mother_is_eta(pi0_has_eta);
      }
      else if (pid == 221)
      {
        bool eta_3pi0 = false;
        bool eta_pi0pipm = false;
        classify_eta_decay(bc, /*hepMCmom*/ nullptr, eta_3pi0, eta_pi0pipm);
        mes->set_eta_3pi0(eta_3pi0);
        mes->set_eta_pi0pipm(eta_pi0pipm);
      }

      std::vector<PHG4Particle *> gamsG4;
      auto itKids = g4_children.find(p->get_track_id());
      if (itKids == g4_children.end())
      {
        continue;
      }

      for (auto *c : itKids->second)
      {
        if (!c)
        {
          continue;
        }
        if (std::abs(c->get_pid()) != 22)
        {
          continue;
        }
        gamsG4.push_back(c);
      }

      if (Verbosity() > 0)
      {
        if (gamsG4.size() > 2)
        {
          std::cout << " In G4 loop more than 2 photon decay check..." << std::endl;
          for (auto *gph : gamsG4)
          {
            std::cout << " pid " << gph->get_pid() << ", parent id " << pid << std::endl;
          }
        }
      }

      if (gamsG4.size() == 2)
      {
        for (auto *gph : gamsG4)
        {
          float ge;
          float gpt;
          float geta;
          float gphi;
          float gpmom;
          fill_kin_from_g4(gph, ge, gpt, geta, gphi, gpmom);
          bool isconverted = false;
          int gphtrackid = gph->get_track_id();
          isconverted = FindConversion(truthinfo, gphtrackid, ge);
          mes->add_photon(ge, gpt, geta, gphi, gpmom, isconverted);
        }
      }
      _container->AddMeson(mes);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void TruthNeutralMesonBuilder::classify_eta_decay(int bc, HepMC::GenParticle *hepMom, bool &eta_3pi0, bool &eta_pi0pipm)
{
  eta_3pi0 = false;
  eta_pi0pipm = false;

  if (hepMom)
  {
    HepMC::GenVertex *vtx = hepMom->end_vertex();
    if (vtx)
    {
      int nPi0 = 0;
      int nPip = 0;
      int nPim = 0;
      for (auto it = vtx->particles_out_const_begin(); it != vtx->particles_out_const_end(); ++it)
      {
        const HepMC::GenParticle *ch = *it;
        if (!ch)
        {
          continue;
        }
        const int chId = ch->pdg_id();
        if (chId == 111)
        {
          ++nPi0;
        }
        else if (chId == 211)
        {
          ++nPip;
        }
        else if (chId == -211)
        {
          ++nPim;
        }
      }
      if (nPi0 == 3)
      {
        eta_3pi0 = true;
      }
      if (nPi0 == 1 && nPip == 1 && nPim == 1)
      {
        eta_pi0pipm = true;
      }
    }
  }

  if (!eta_3pi0 && !eta_pi0pipm && truthinfo)
  {
    auto itMom = g4_by_barcode.find(bc);
    if (itMom != g4_by_barcode.end())
    {
      PHG4Particle *g4mom = itMom->second;
      auto itKids = g4_children.find(g4mom->get_track_id());
      if (itKids != g4_children.end())
      {
        int nPi0 = 0;
        int nPip = 0;
        int nPim = 0;
        for (auto *c : itKids->second)
        {
          if (!c)
          {
            continue;
          }
          const int chId = c->get_pid();
          if (chId == 111)
          {
            ++nPi0;
          }
          else if (chId == 211)
          {
            ++nPip;
          }
          else if (chId == -211)
          {
            ++nPim;
          }
        }
        if (nPi0 == 3)
        {
          eta_3pi0 = true;
        }
        if (nPi0 == 1 && nPip == 1 && nPim == 1)
        {
          eta_pi0pipm = true;
        }
      }
    }
  }
}

bool TruthNeutralMesonBuilder::FindConversion(PHG4TruthInfoContainer *_truthinfo, int trackid, float energy)
{
  PHG4Shower *shower = _truthinfo->GetShower(trackid);
  if (!shower)
  {
    return false;
  }
  bool foundconversion = false;
  auto g4particle_ids = shower->g4particle_ids();
  for (auto g4particle_id : g4particle_ids)
  {
    PHG4Particle *g4particle = _truthinfo->GetParticle(g4particle_id);
    if (!g4particle)
    {
      continue;
    }
    int vertexid = g4particle->get_vtx_id();
    PHG4VtxPoint *vtxp = _truthinfo->GetVtx(vertexid);
    if (!vtxp)
    {
      continue;
    }
    float vertexr = sqrt(vtxp->get_x() * vtxp->get_x() + vtxp->get_y() * vtxp->get_y());
    if (vertexr < m_conversion_radius_limit)
    {
      float momentum = sqrt(g4particle->get_px() * g4particle->get_px() + g4particle->get_py() * g4particle->get_py() + g4particle->get_pz() * g4particle->get_pz());
      if (momentum > 0.3 * energy)
      {
        int g4particlepid = g4particle->get_pid();
        if (abs(g4particlepid) == 11)
        {
          foundconversion = true;
        }
      }
    }
  }
  return foundconversion;
}

bool TruthNeutralMesonBuilder::RejectShowerMeson(PHG4TruthInfoContainer *_truthinfo, int parent_trackid, int this_trackid)
{
  PHG4Shower *shower = _truthinfo->GetShower(parent_trackid);
  if (!shower)
  {
    return false;
  }
  bool reject = false;
  auto g4particle_ids = shower->g4particle_ids();
  for (auto g4particle_id : g4particle_ids)
  {
    if (g4particle_id != this_trackid)
    {
      continue;
    }
    PHG4Particle *g4particle = _truthinfo->GetParticle(g4particle_id);
    if (!g4particle)
    {
      continue;
    }
    int vertexid = g4particle->get_vtx_id();
    PHG4VtxPoint *vtxp = _truthinfo->GetVtx(vertexid);
    if (!vtxp)
    {
      continue;
    }

    float vertexr = sqrt(vtxp->get_x() * vtxp->get_x() + vtxp->get_y() * vtxp->get_y());
    if (vertexr > m_shower_reject_radius)
    {
      int g4particlepid = g4particle->get_pid();
      if (abs(g4particlepid) == 111 || abs(g4particlepid) == 221)
      {
        reject = true;
      }
    }
  }
  return reject;
}

int TruthNeutralMesonBuilder::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
