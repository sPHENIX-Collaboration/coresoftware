#include "SvtxTrack_v5.h"
#include "SvtxTrackState.h"
#include "SvtxTrackState_v1.h"
#include "TrackSeed.h"

#include <trackbase/TrkrDefs.h>

#include <phool/PHObject.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <utility>

namespace
{
  void count_cluster_key(TrkrDefs::cluskey key,
                         unsigned int& nmvtx,
                         unsigned int& nintt,
                         unsigned int& ntpc,
                         unsigned int& ntpot)
  {
    switch (TrkrDefs::getTrkrId(key))
    {
    case TrkrDefs::mvtxId:
      ++nmvtx;
      break;
    case TrkrDefs::inttId:
      ++nintt;
      break;
    case TrkrDefs::tpcId:
      ++ntpc;
      break;
    case TrkrDefs::micromegasId:
      ++ntpot;
      break;
    default:
      break;
    }
  }

  void count_seed_clusters(const TrackSeed* seed,
                           unsigned int& nmvtx,
                           unsigned int& nintt,
                           unsigned int& ntpc,
                           unsigned int& ntpot)
  {
    if (!seed)
    {
      return;
    }

    for (auto iter = seed->begin_cluster_keys(); iter != seed->end_cluster_keys(); ++iter)
    {
      count_cluster_key(*iter, nmvtx, nintt, ntpc, ntpot);
    }
  }
}

SvtxTrack_v5::SvtxTrack_v5()
{
  // always include the pca point
  _states.insert(std::make_pair(0, new SvtxTrackState_v1(0)));
}

SvtxTrack_v5::SvtxTrack_v5(const SvtxTrack& source)
{
  SvtxTrack_v5::CopyFrom(source);
}

// have to suppress missingMemberCopy from cppcheck, it does not
// go down to the CopyFrom method where things are done correctly
// cppcheck-suppress missingMemberCopy
SvtxTrack_v5::SvtxTrack_v5(const SvtxTrack_v5& source)
  : SvtxTrack(source)
{
  SvtxTrack_v5::CopyFrom(source);
}

SvtxTrack_v5& SvtxTrack_v5::operator=(const SvtxTrack_v5& source)
{
  if (this != &source)
  {
    CopyFrom(source);
  }
  return *this;
}

SvtxTrack_v5::~SvtxTrack_v5()
{
  clear_states();
}

void SvtxTrack_v5::CopyFrom(const SvtxTrack& source)
{
  // do nothing if copying onto oneself
  if (this == &source)
  {
    return;
  }

  // parent class method
  SvtxTrack::CopyFrom(source);

  _track_id = source.get_id();
  _vertex_id = source.get_vertex_id();
  m_charge = static_cast<signed char>((source.get_charge() > 0) ? 1 : -1);
  _chisq = source.get_chisq();
  set_ndf(source.get_ndf());
  _track_crossing = source.get_crossing();

  clear_states();
  for (auto iter = source.begin_states(); iter != source.end_states(); ++iter)
  {
    _states.insert(std::make_pair(iter->first, static_cast<SvtxTrackState*>(iter->second->CloneMe())));
  }

  const auto* source_v5 = dynamic_cast<const SvtxTrack_v5*>(&source);
  if (source_v5)
  {
    set_nmvtx_clusters(source_v5->get_nmvtx_clusters());
    set_nintt_clusters(source_v5->get_nintt_clusters());
    set_ntpc_clusters(source_v5->get_ntpc_clusters());
    set_ntpot_clusters(source_v5->get_ntpot_clusters());
    return;
  }

  unsigned int nmvtx = 0;
  unsigned int nintt = 0;
  unsigned int ntpc = 0;
  unsigned int ntpot = 0;

  const auto begin_cluster_keys = source.begin_cluster_keys();
  const auto end_cluster_keys = source.end_cluster_keys();
  if (begin_cluster_keys != end_cluster_keys)
  {
    for (auto iter = begin_cluster_keys; iter != end_cluster_keys; ++iter)
    {
      count_cluster_key(*iter, nmvtx, nintt, ntpc, ntpot);
    }
  }
  else
  {
    count_seed_clusters(source.get_silicon_seed(), nmvtx, nintt, ntpc, ntpot);
    count_seed_clusters(source.get_tpc_seed(), nmvtx, nintt, ntpc, ntpot);
  }

  set_nmvtx_clusters(nmvtx);
  set_nintt_clusters(nintt);
  set_ntpc_clusters(ntpc);
  set_ntpot_clusters(ntpot);
}

void SvtxTrack_v5::identify(std::ostream& os) const
{
  os << "SvtxTrack_v5 Object ";
  os << "id: " << get_id() << " ";
  os << "vertex id: " << get_vertex_id() << " ";
  os << "charge: " << get_charge() << " ";
  os << "chisq: " << get_chisq() << " ndf:" << get_ndf() << " ";
  os << "crossing: " << get_crossing() << " ";
  os << "clusters: MVTX " << get_nmvtx_clusters()
     << " INTT " << get_nintt_clusters()
     << " TPC " << get_ntpc_clusters()
     << " TPOT " << get_ntpot_clusters() << " ";
  os << "nstates: " << _states.size() << " ";
  os << std::endl;

  os << "(px,py,pz) = ("
     << get_px() << ","
     << get_py() << ","
     << get_pz() << ")" << std::endl;

  os << "(x,y,z) = (" << get_x() << "," << get_y() << "," << get_z() << ")" << std::endl;
}

int SvtxTrack_v5::isValid() const
{
  return 1;
}

void SvtxTrack_v5::set_ndf(int ndf)
{
  _ndf = static_cast<unsigned char>(std::max(0, std::min(ndf, static_cast<int>(std::numeric_limits<unsigned char>::max()))));
}

float SvtxTrack_v5::get_x() const
{
  const auto* state = get_pca_state();
  return state ? state->get_x() : NAN;
}

void SvtxTrack_v5::set_x(float x)
{
  get_pca_state()->set_x(x);
}

float SvtxTrack_v5::get_y() const
{
  const auto* state = get_pca_state();
  return state ? state->get_y() : NAN;
}

void SvtxTrack_v5::set_y(float y)
{
  get_pca_state()->set_y(y);
}

float SvtxTrack_v5::get_z() const
{
  const auto* state = get_pca_state();
  return state ? state->get_z() : NAN;
}

void SvtxTrack_v5::set_z(float z)
{
  get_pca_state()->set_z(z);
}

float SvtxTrack_v5::get_pos(unsigned int i) const
{
  const auto* state = get_pca_state();
  return state ? state->get_pos(i) : NAN;
}

float SvtxTrack_v5::get_px() const
{
  const auto* state = get_pca_state();
  return state ? state->get_px() : NAN;
}

void SvtxTrack_v5::set_px(float px)
{
  get_pca_state()->set_px(px);
}

float SvtxTrack_v5::get_py() const
{
  const auto* state = get_pca_state();
  return state ? state->get_py() : NAN;
}

void SvtxTrack_v5::set_py(float py)
{
  get_pca_state()->set_py(py);
}

float SvtxTrack_v5::get_pz() const
{
  const auto* state = get_pca_state();
  return state ? state->get_pz() : NAN;
}

void SvtxTrack_v5::set_pz(float pz)
{
  get_pca_state()->set_pz(pz);
}

float SvtxTrack_v5::get_mom(unsigned int i) const
{
  const auto* state = get_pca_state();
  return state ? state->get_mom(i) : NAN;
}

float SvtxTrack_v5::get_error(int i, int j) const
{
  const auto* state = get_pca_state();
  return state ? state->get_error(i, j) : NAN;
}

void SvtxTrack_v5::set_error(int i, int j, float value)
{
  get_pca_state()->set_error(i, j, value);
}

void SvtxTrack_v5::clear_states()
{
  for (const auto& pair : _states)
  {
    delete pair.second;
  }

  _states.clear();
}

const SvtxTrackState* SvtxTrack_v5::get_state(float pathlength) const
{
  const auto iter = _states.find(pathlength);
  return (iter == _states.end()) ? nullptr : iter->second;
}

SvtxTrackState* SvtxTrack_v5::get_state(float pathlength)
{
  const auto iter = _states.find(pathlength);
  return (iter == _states.end()) ? nullptr : iter->second;
}

SvtxTrackState* SvtxTrack_v5::insert_state(const SvtxTrackState* state)
{
  if (!state)
  {
    return nullptr;
  }

  const auto pathlength = state->get_pathlength();
  auto iterator = _states.lower_bound(pathlength);
  if (iterator == _states.end() || pathlength < iterator->first)
  {
    auto* const copy = static_cast<SvtxTrackState*>(state->CloneMe());
    iterator = _states.insert(iterator, std::make_pair(pathlength, copy));
  }

  return iterator->second;
}

size_t SvtxTrack_v5::erase_state(float pathlength)
{
  StateIter iter = _states.find(pathlength);
  if (iter == _states.end())
  {
    return _states.size();
  }

  delete iter->second;
  _states.erase(iter);
  return _states.size();
}

unsigned char SvtxTrack_v5::compress_cluster_count(unsigned int nclusters)
{
  return static_cast<unsigned char>(std::min(nclusters, static_cast<unsigned int>(std::numeric_limits<unsigned char>::max())));
}

const SvtxTrackState* SvtxTrack_v5::get_pca_state() const
{
  const auto iter = _states.find(0.0);
  return (iter == _states.end()) ? nullptr : iter->second;
}

SvtxTrackState* SvtxTrack_v5::get_pca_state()
{
  auto iter = _states.find(0.0);
  if (iter == _states.end())
  {
    iter = _states.insert(std::make_pair(0.0, new SvtxTrackState_v1(0.0))).first;
  }

  return iter->second;
}
