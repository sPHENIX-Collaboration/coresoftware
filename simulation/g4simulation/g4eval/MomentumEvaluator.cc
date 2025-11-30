#include "MomentumEvaluator.h"

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TSystem.h>

#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <utility>
#include <vector>

class TrivialTrack
{
 public:
  // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
  float px, py, pz;
  float dcax, dcay, dcaz;
  float quality;
  // NOLINTEND(misc-non-private-member-variables-in-classes)
  TrivialTrack(float x, float y, float z, float dx, float dy, float dz, float qual = 0.)
    : px(x)
    , py(y)
    , pz(z)
    , dcax(dx)
    , dcay(dy)
    , dcaz(dz)
    , quality(qual)
  {
  }
};

class RecursiveMomentumContainer  // NOLINT(hicpp-special-member-functions)
{
 public:
  RecursiveMomentumContainer(float PX_LO, float PX_HI, float PY_LO, float PY_HI, float PZ_LO, float PZ_HI, int MLEV, int LEV = 0)
    : px_lo(PX_LO)
    , px_hi(PX_HI)
    , py_lo(PY_LO)
    , py_hi(PY_HI)
    , pz_lo(PZ_LO)
    , pz_hi(PZ_HI)
    , level(LEV)
    , maxlevel(MLEV)
    , x_pos(0)
    , y_pos(0)
    , z_pos(0)
  {
    for (auto& container : containers)
    {
      for (auto& j : container)
      {
        for (auto& k : j)
        {
          k = nullptr;
        }
      }
    }
  }
  virtual ~RecursiveMomentumContainer()
  {
    for (auto& container : containers)
    {
      for (auto& j : container)
      {
        for (auto& k : j)
        {
          delete k;
        }
      }
    }
  }

  virtual bool insert(TrivialTrack& track);

  virtual TrivialTrack* begin()
  {
    x_pos = 0;
    y_pos = 0;
    z_pos = 0;
    while (true)
    {
      if (containers[x_pos][y_pos][z_pos] == nullptr)
      {
        if (z_pos == 0)
        {
          z_pos = 1;
          continue;
        }

        if (y_pos == 0)
        {
          z_pos = 0;
          y_pos = 1;
          continue;
        }

        if (x_pos == 0)
        {
          z_pos = 0;
          y_pos = 0;
          x_pos = 1;
          continue;
        }
        return nullptr;
      }
      return containers[x_pos][y_pos][z_pos]->begin();
    }
  }

  virtual TrivialTrack* next()
  {
    bool block_changed = false;
    while (true)
    {
      if (containers[x_pos][y_pos][z_pos] == nullptr)
      {
        block_changed = true;
        if (z_pos == 0)
        {
          z_pos = 1;
          continue;
        }

        if (y_pos == 0)
        {
          z_pos = 0;
          y_pos = 1;
          continue;
        }

        if (x_pos == 0)
        {
          z_pos = 0;
          y_pos = 0;
          x_pos = 1;
          continue;
        }

        return nullptr;
      }
      TrivialTrack* val = nullptr;
      if (block_changed == true)
      {
        val = containers[x_pos][y_pos][z_pos]->begin();
      }
      else
      {
        val = containers[x_pos][y_pos][z_pos]->next();
      }

      if (val == nullptr)
      {
        block_changed = true;
        if (z_pos == 0)
        {
          z_pos = 1;
          continue;
        }

        if (y_pos == 0)
        {
          z_pos = 0;
          y_pos = 1;
          continue;
        }

        if (x_pos == 0)
        {
          z_pos = 0;
          y_pos = 0;
          x_pos = 1;
          continue;
        }
        return nullptr;
      }

      return val;
    }
  }

  virtual void append_list(std::vector<TrivialTrack*>& track_list, float PX_LO, float PX_HI, float PY_LO, float PY_HI, float PZ_LO, float PZ_HI)
  {
    for (auto& container : containers)
    {
      for (auto& j : container)
      {
        for (auto& k : j)
        {
          if (k == nullptr)
          {
            continue;
          }

          if ((k->px_hi < PX_LO) || (k->px_lo > PX_HI) || (k->py_hi < PY_LO) || (k->py_lo > PY_HI) || (k->pz_hi < PZ_LO) || (k->pz_lo > PZ_HI))
          {
            continue;
          }

          k->append_list(track_list, PX_LO, PX_HI, PY_LO, PY_HI, PZ_LO, PZ_HI);
        }
      }
    }
  }

 protected:
  // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
  float px_lo, px_hi;
  float py_lo, py_hi;
  float pz_lo, pz_hi;

  int level;
  int maxlevel;

  unsigned int x_pos, y_pos, z_pos;

  RecursiveMomentumContainer* containers[2][2][2]{};
  // NOLINTEND(misc-non-private-member-variables-in-classes)
};

class RecursiveMomentumContainerEnd : public RecursiveMomentumContainer
{
 public:
  RecursiveMomentumContainerEnd(float PX_LO, float PX_HI, float PY_LO, float PY_HI, float PZ_LO, float PZ_HI, int MLEV, int LEV = 0)
    : RecursiveMomentumContainer(PX_LO, PX_HI, PY_LO, PY_HI, PZ_LO, PZ_HI, MLEV, LEV)
  {
  }

  //  ~RecursiveMomentumContainerEnd() override = default;

  bool insert(TrivialTrack& track) override
  {
    tracks.push_back(track);
    return true;
  }

  TrivialTrack* begin() override
  {
    x_pos = 0;
    return (&(tracks.at(0)));
  }

  TrivialTrack* next() override
  {
    if (x_pos >= (tracks.size() - 1))
    {
      return nullptr;
    }

    x_pos += 1;
    return (&(tracks[x_pos]));
  }

  void append_list(std::vector<TrivialTrack*>& track_list, float PX_LO, float PX_HI, float PY_LO, float PY_HI, float PZ_LO, float PZ_HI) override
  {
    for (auto& track : tracks)
    {
      if ((track.px < PX_LO) || (track.px > PX_HI) || (track.py < PY_LO) || (track.py > PY_HI) || (track.pz < PZ_LO) || (track.pz > PZ_HI))
      {
        continue;
      }
      track_list.push_back(&track);
    }
  }

 protected:
  std::vector<TrivialTrack> tracks;  // NOLINT(misc-non-private-member-variables-in-classes)
};

bool RecursiveMomentumContainer::insert(TrivialTrack& track)
{
  if ((track.px < px_lo) || (track.py < py_lo) || (track.pz < pz_lo) || (track.px > px_hi) || (track.py > py_hi) || (track.pz > pz_hi))
  {
    return false;
  }

  int x_ind = 0;
  if (track.px > (px_lo + 0.5 * (px_hi - px_lo)))
  {
    x_ind = 1;
  }
  int y_ind = 0;
  if (track.py > (py_lo + 0.5 * (py_hi - py_lo)))
  {
    y_ind = 1;
  }
  int z_ind = 0;
  if (track.pz > (pz_lo + 0.5 * (pz_hi - pz_lo)))
  {
    z_ind = 1;
  }

  if (containers[x_ind][y_ind][z_ind] == nullptr)
  {
    float px_lo_new = px_lo + (float(x_ind)) * 0.5 * (px_hi - px_lo);
    float px_hi_new = px_lo_new + 0.5 * (px_hi - px_lo);

    float py_lo_new = py_lo + (float(y_ind)) * 0.5 * (py_hi - py_lo);
    float py_hi_new = py_lo_new + 0.5 * (py_hi - py_lo);

    float pz_lo_new = pz_lo + (float(z_ind)) * 0.5 * (pz_hi - pz_lo);
    float pz_hi_new = pz_lo_new + 0.5 * (pz_hi - pz_lo);

    if (level < maxlevel)
    {
      containers[x_ind][y_ind][z_ind] = new RecursiveMomentumContainer(px_lo_new, px_hi_new, py_lo_new, py_hi_new, pz_lo_new, pz_hi_new, maxlevel, level + 1);
    }
    else
    {
      containers[x_ind][y_ind][z_ind] = new RecursiveMomentumContainerEnd(px_lo_new, px_hi_new, py_lo_new, py_hi_new, pz_lo_new, pz_hi_new, maxlevel, level + 1);
    }
  }
  return containers[x_ind][y_ind][z_ind]->insert(track);
}

MomentumEvaluator::MomentumEvaluator(const std::string& fname, float pt_s, float pz_s, unsigned int /*n_l*/, unsigned int n_i, unsigned int n_r, float i_z, float o_z)
  : ntp_true(nullptr)
  , ntp_reco(nullptr)
  , pt_search_scale(pt_s)
  , pz_search_scale(pz_s)
  , event_counter(0)
  , file_name(fname)
  , n_inner_layers(n_i)
  , n_required_layers(n_r)
  , inner_z_length(i_z)
  , outer_z_length(o_z)
{
}
MomentumEvaluator::~MomentumEvaluator()
{
  delete ntp_true;

  delete ntp_reco;
}

int MomentumEvaluator::Init(PHCompositeNode* /*topNode*/)
{
  delete ntp_true;

  delete ntp_reco;

  ntp_true = new TNtuple("ntp_true", "true simulated tracks", "event:px:py:pz:dcax:dcay:dcaz:r_px:r_py:r_pz:r_dcax:r_dcay:r_dcaz:quality");
  ntp_reco = new TNtuple("ntp_reco", "reconstructed tracks", "event:px:py:pz:dcax:dcay:dcaz:t_px:t_py:t_pz:t_dcax:t_dcay:t_dcaz:quality");
  event_counter = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

int MomentumEvaluator::process_event(PHCompositeNode* topNode)
{
  PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  PHG4HitContainer* g4hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_SVTX");
  if (g4hits == nullptr)
  {
    std::cout << "can't find PHG4HitContainer" << std::endl;
    exit(1);
  }
  PHG4HitContainer::ConstRange g4range = g4hits->getHits();

  // set<int> trkids;
  std::map<int, std::pair<unsigned int, unsigned int> > trkids;

  for (PHG4HitContainer::ConstIterator iter = g4range.first; iter != g4range.second; ++iter)
  {
    PHG4Hit* hit = iter->second;

    int layer = hit->get_layer();
    float length = outer_z_length;
    if (((unsigned int) layer) < n_inner_layers)
    {
      length = inner_z_length;
    }
    if (std::abs(hit->get_z(0)) > length)
    {
      continue;
    }

    int trk_id = hit->get_trkid();
    if (!trkids.contains(trk_id))
    {
      trkids[trk_id].first = 0;
      trkids[trk_id].second = 0;
    }
    if (hit->get_layer() < 32)
    {
      trkids[trk_id].first = (trkids[trk_id].first | (1U << (hit->get_layer())));
    }
    else if (hit->get_layer() < 64)
    {
      trkids[trk_id].second = (trkids[trk_id].second | (1U << (hit->get_layer() - 32U)));
    }
    else
    {
      std::cout << PHWHERE << "hit layer out of bounds (0-63) " << hit->get_layer() << std::endl;
      gSystem->Exit(1);
      exit(1);
    }

    // std::cout<<"trk_id = "<<trk_id<<std::endl;
    // std::cout<<"layer = "<<hit->get_layer()<<std::endl;
    // std::cout<<"nlayer = "<<__builtin_popcount(trkids[trk_id].first)+__builtin_popcount(trkids[trk_id].second)<<std::endl<<std::endl;
    // trkids.insert(trk_id);
  }

  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  PHG4VtxPoint* gvertex = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
  float gvx = gvertex->get_x();
  float gvy = gvertex->get_y();
  float gvz = gvertex->get_z();

  RecursiveMomentumContainer true_sorted(-20., 20., -20., 20., -20., 20., 10);

  // PHG4TruthInfoContainer::Map primarymap = truthinfo->GetPrimaryMap();
  PHG4TruthInfoContainer::Map primarymap = truthinfo->GetMap();
  for (auto& iter : primarymap)
  {
    PHG4Particle* particle = iter.second;

    float vx = truthinfo->GetVtx(particle->get_vtx_id())->get_x();
    float vy = truthinfo->GetVtx(particle->get_vtx_id())->get_y();
    float vz = truthinfo->GetVtx(particle->get_vtx_id())->get_z();

    TrivialTrack track(particle->get_px(), particle->get_py(), particle->get_pz(), vx - gvx, vy - gvy, vz - gvz);

    if (((track.px * track.px) + (track.py * track.py)) < (0.1 * 0.1))
    {
      continue;
    }

    if (!trkids.contains(particle->get_track_id()))
    {
      continue;
    }

    // std::cout<<"trk, nhits = "<<particle->get_track_id()<<" "<<__builtin_popcount(trkids[particle->get_track_id()].first)+__builtin_popcount(trkids[particle->get_track_id()].second)<<std::endl;

    if (__builtin_popcount(trkids[particle->get_track_id()].first) + __builtin_popcount(trkids[particle->get_track_id()].second) < (int) n_required_layers)
    {
      continue;
    }

    true_sorted.insert(track);
  }

  RecursiveMomentumContainer reco_sorted(-20., 20., -20., 20., -20., 20., 10);
  for (auto& iter : *trackmap)
  {
    SvtxTrack* track = iter.second;

    TrivialTrack ttrack(track->get_px(), track->get_py(), track->get_pz(), track->get_x() - gvx, track->get_y() - gvy, track->get_z() - gvz, track->get_quality());
    reco_sorted.insert(ttrack);
  }

  TrivialTrack* t_track = true_sorted.begin();
  std::vector<TrivialTrack*> pointer_list;
  while (t_track != nullptr)
  {
    pointer_list.clear();

    float pt = std::sqrt((t_track->px * t_track->px) + (t_track->py * t_track->py));
    float pt_diff = pt * pt_search_scale;
    float px_lo = t_track->px - pt_diff;
    float px_hi = t_track->px + pt_diff;
    float py_lo = t_track->py - pt_diff;
    float py_hi = t_track->py + pt_diff;
    float pz_diff = std::abs(t_track->pz) * pz_search_scale;
    float pz_lo = t_track->pz - pz_diff;
    float pz_hi = t_track->pz + pz_diff;

    reco_sorted.append_list(pointer_list, px_lo, px_hi, py_lo, py_hi, pz_lo, pz_hi);

    if (!pointer_list.empty())
    {
      float mom_true = std::sqrt(pt * pt + (t_track->pz) * (t_track->pz));
      float best_ind = 0;
      float mom_reco = std::sqrt((pointer_list[0]->px) * (pointer_list[0]->px) + (pointer_list[0]->py) * (pointer_list[0]->py) + (pointer_list[0]->pz) * (pointer_list[0]->pz));
      float best_mom = mom_reco;
      for (unsigned int i = 1; i < pointer_list.size(); ++i)
      {
        mom_reco = std::sqrt((pointer_list[i]->px) * (pointer_list[i]->px) + (pointer_list[i]->py) * (pointer_list[i]->py) + (pointer_list[i]->pz) * (pointer_list[i]->pz));
        if (std::abs(mom_true - mom_reco) < std::abs(mom_true - best_mom))
        {
          best_mom = mom_reco;
          best_ind = i;
        }
      }

      float ntp_data[14] = {(float) event_counter, t_track->px, t_track->py, t_track->pz, t_track->dcax, t_track->dcay, t_track->dcaz, pointer_list[best_ind]->px, pointer_list[best_ind]->py, pointer_list[best_ind]->pz, pointer_list[best_ind]->dcax, pointer_list[best_ind]->dcay, pointer_list[best_ind]->dcaz, pointer_list[best_ind]->quality};
      ntp_true->Fill(ntp_data);
    }
    else
    {
      float ntp_data[14] = {(float) event_counter, t_track->px, t_track->py, t_track->pz, t_track->dcax, t_track->dcay, t_track->dcaz, -9999., -9999., -9999., -9999., -9999., -9999., -9999.};
      ntp_true->Fill(ntp_data);
    }

    t_track = true_sorted.next();
  }

  TrivialTrack* r_track = reco_sorted.begin();
  while (r_track != nullptr)
  {
    pointer_list.clear();

    float pt = std::sqrt((r_track->px * r_track->px) + (r_track->py * r_track->py));
    float pt_diff = pt * pt_search_scale;
    float px_lo = r_track->px - pt_diff;
    float px_hi = r_track->px + pt_diff;
    float py_lo = r_track->py - pt_diff;
    float py_hi = r_track->py + pt_diff;
    float pz_diff = std::abs(r_track->pz) * pz_search_scale;
    float pz_lo = r_track->pz - pz_diff;
    float pz_hi = r_track->pz + pz_diff;

    true_sorted.append_list(pointer_list, px_lo, px_hi, py_lo, py_hi, pz_lo, pz_hi);

    if (!pointer_list.empty())
    {
      float mom_reco = std::sqrt(pt * pt + (r_track->pz) * (r_track->pz));
      float best_ind = 0;
      float mom_true = std::sqrt((pointer_list[0]->px) * (pointer_list[0]->px) + (pointer_list[0]->py) * (pointer_list[0]->py) + (pointer_list[0]->pz) * (pointer_list[0]->pz));
      float best_mom = mom_true;
      for (unsigned int i = 1; i < pointer_list.size(); ++i)
      {
        mom_true = std::sqrt((pointer_list[i]->px) * (pointer_list[i]->px) + (pointer_list[i]->py) * (pointer_list[i]->py) + (pointer_list[i]->pz) * (pointer_list[i]->pz));
        if (std::abs(mom_reco - mom_true) < std::abs(mom_reco - best_mom))
        {
          best_mom = mom_true;
          best_ind = i;
        }
      }

      float ntp_data[14] = {(float) event_counter, r_track->px, r_track->py, r_track->pz, r_track->dcax, r_track->dcay, r_track->dcaz, pointer_list[best_ind]->px, pointer_list[best_ind]->py, pointer_list[best_ind]->pz, pointer_list[best_ind]->dcax, pointer_list[best_ind]->dcay, pointer_list[best_ind]->dcaz, r_track->quality};
      ntp_reco->Fill(ntp_data);
    }
    else
    {
      float ntp_data[14] = {(float) event_counter, r_track->px, r_track->py, r_track->pz, r_track->dcax, r_track->dcay, r_track->dcaz, -9999., -9999., -9999., -9999., -9999., -9999., r_track->quality};
      ntp_reco->Fill(ntp_data);
    }

    r_track = reco_sorted.next();
  }

  event_counter += 1;
  return Fun4AllReturnCodes::EVENT_OK;
}

int MomentumEvaluator::End(PHCompositeNode* /*topNode*/)
{
  TFile outfile(file_name.c_str(), "recreate");
  outfile.cd();
  ntp_true->Write();
  ntp_reco->Write();
  outfile.Close();

  return Fun4AllReturnCodes::EVENT_OK;
}
