#include "SvtxTrack_v1.h"
#include "SvtxTrackState.h"
#include "SvtxTrackState_v1.h"

#include <trackbase/TrkrDefs.h>  // for cluskey

#include <cassert>
#include <climits>
#include <map>
#include <vector>                // for vector


using namespace std;

SvtxTrack_v1::SvtxTrack_v1()
  : _track_id(UINT_MAX)
  , _vertex_id(UINT_MAX)
  , _is_positive_charge(false)
  , _chisq(NAN)
  , _ndf(0)
  , _dca(NAN)
  , _dca_error(NAN)
  , _dca2d(NAN)
  , _dca2d_error(NAN)
  , _dca3d_xy(NAN)
  , _dca3d_xy_error(NAN)
  , _dca3d_z(NAN)
  , _dca3d_z_error(NAN)
  , _states()
  , _cluster_ids()
  , _cluster_keys()
  , _cal_dphi()
  , _cal_deta()
  , _cal_energy_3x3()
  , _cal_energy_5x5()
  , _cal_cluster_id()
  , _cal_cluster_key()
  , _cal_cluster_e()
{
  // always include the pca point
  _states.insert(make_pair(0.0, new SvtxTrackState_v1(0.0)));
}

SvtxTrack_v1::SvtxTrack_v1(const SvtxTrack_v1& track)
{
  *this = track;
  return;
}

SvtxTrack_v1& SvtxTrack_v1::operator=(const SvtxTrack_v1& track)
{
  _track_id = track.get_id();
  _vertex_id = track.get_vertex_id();
  _is_positive_charge = track.get_positive_charge();
  _chisq = track.get_chisq();
  _ndf = track.get_ndf();
  _dca = track.get_dca();
  _dca_error = track.get_dca_error();
  _dca2d = track.get_dca2d();
  _dca2d_error = track.get_dca2d_error();
  _dca3d_xy = track.get_dca3d_xy();
  _dca3d_xy_error = track.get_dca3d_xy_error();
  _dca3d_z = track.get_dca3d_z();
  _dca3d_z_error = track.get_dca3d_z_error();

  // copy the states over into new state objects stored here
  clear_states();
  for (ConstStateIter iter = track.begin_states();
       iter != track.end_states();
       ++iter)
  {
    SvtxTrackState* state = dynamic_cast< SvtxTrackState*> (iter->second->CloneMe());
    _states.insert(make_pair(state->get_pathlength(), state));
  }

  // copy over cluster ID set
  _cluster_ids.clear();
  for (ConstClusterIter iter = track.begin_clusters();
       iter != track.end_clusters();
       ++iter)
  {
    _cluster_ids.insert(*iter);
  }

  // copy over cluster key set
  _cluster_keys.clear();
  for (ConstClusterKeyIter iter = track.begin_cluster_keys();
       iter != track.end_cluster_keys();
       ++iter)
  {
    _cluster_keys.insert(*iter);
  }

  // copy over calorimeter projections
  std::vector<CAL_LAYER> types;
  types.push_back(SvtxTrack::PRES);
  types.push_back(SvtxTrack::CEMC);
  types.push_back(SvtxTrack::HCALIN);
  types.push_back(SvtxTrack::HCALOUT);

  _cal_dphi.clear();
  _cal_deta.clear();
  _cal_energy_3x3.clear();
  _cal_energy_5x5.clear();
  _cal_cluster_id.clear();
  _cal_cluster_key.clear();
  _cal_cluster_e.clear();

  for (unsigned int i = 0; i < types.size(); ++i)
  {
    if (!isnan(track.get_cal_dphi(types[i]))) set_cal_dphi(types[i], track.get_cal_dphi(types[i]));
    if (!isnan(track.get_cal_deta(types[i]))) set_cal_deta(types[i], track.get_cal_deta(types[i]));
    if (!isnan(track.get_cal_energy_3x3(types[i]))) set_cal_energy_3x3(types[i], track.get_cal_energy_3x3(types[i]));
    if (!isnan(track.get_cal_energy_5x5(types[i]))) set_cal_energy_5x5(types[i], track.get_cal_energy_5x5(types[i]));
    if (track.get_cal_cluster_id(types[i]) != UINT_MAX) set_cal_cluster_id(types[i], track.get_cal_cluster_id(types[i]));
    if (track.get_cal_cluster_key(types[i]) != UINT_MAX) set_cal_cluster_key(types[i], track.get_cal_cluster_key(types[i]));
    if (!isnan(track.get_cal_cluster_e(types[i]))) set_cal_cluster_e(types[i], track.get_cal_cluster_e(types[i]));
  }

  return *this;
}

SvtxTrack_v1::~SvtxTrack_v1()
{
  clear_states();
}

void SvtxTrack_v1::identify(std::ostream& os) const
{
  os << "SvtxTrack_v1 Object ";
  os << "id: " << get_id() << " ";
  os << "vertex id: " << get_vertex_id() << " ";
  os << "charge: " << get_charge() << " ";
  os << "chisq: " << get_chisq() << " ndf:" << get_ndf() << " ";
  os << endl;

  os << "(px,py,pz) = ("
     << get_px() << ","
     << get_py() << ","
     << get_pz() << ")" << endl;

  os << "(x,y,z) = (" << get_x() << "," << get_y() << "," << get_z() << ")" << endl;

  if ( _cluster_ids.size() > 0 || _cluster_keys.size() > 0 )
  {
    os << "list of cluster IDs ";
    for (SvtxTrack::ConstClusterIter iter = begin_clusters();
         iter != end_clusters();
         ++iter)
    {
      unsigned int cluster_id = *iter;
      os << cluster_id << " ";
    }

    os << "list of cluster keys ";
    for (SvtxTrack::ConstClusterKeyIter iter = begin_cluster_keys();
         iter != end_cluster_keys();
         ++iter)
    {
      TrkrDefs::cluskey cluster_key = *iter;
      os << cluster_key << " ";
    }
  }
  else
    os << " track has no clusters " << endl;
  
  os << endl;

  return;
}

void SvtxTrack_v1::clear_states()
{
  for (StateIter iter = _states.begin();
       iter != _states.end();
       ++iter)
  {
    SvtxTrackState* state = iter->second;
    delete state;
  }
  _states.clear();
}

int SvtxTrack_v1::isValid() const
{
  return 1;
}

const SvtxTrackState* SvtxTrack_v1::get_state(float pathlength) const
{
  ConstStateIter iter = _states.find(pathlength);
  if (iter == _states.end()) return nullptr;
  return iter->second;
}

SvtxTrackState* SvtxTrack_v1::get_state(float pathlength)
{
  StateIter iter = _states.find(pathlength);
  if (iter == _states.end()) return nullptr;
  return iter->second;
}

SvtxTrackState* SvtxTrack_v1::insert_state(const SvtxTrackState* state)
{
  _states.insert(make_pair(state->get_pathlength(), dynamic_cast< SvtxTrackState*> (state->CloneMe())));
  return _states[state->get_pathlength()];
}

size_t SvtxTrack_v1::erase_state(float pathlength)
{
  StateIter iter = _states.find(pathlength);
  if (iter == _states.end()) return _states.size();

  delete iter->second;
  _states.erase(iter);
  return _states.size();
}

float SvtxTrack_v1::get_cal_dphi(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, float>::const_iterator citer = _cal_dphi.find(layer);
  if (citer == _cal_dphi.end()) return NAN;
  return citer->second;
}

float SvtxTrack_v1::get_cal_deta(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, float>::const_iterator citer = _cal_deta.find(layer);
  if (citer == _cal_deta.end()) return NAN;
  return citer->second;
}

float SvtxTrack_v1::get_cal_energy_3x3(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, float>::const_iterator citer = _cal_energy_3x3.find(layer);
  if (citer == _cal_energy_3x3.end()) return NAN;
  return citer->second;
}

float SvtxTrack_v1::get_cal_energy_5x5(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, float>::const_iterator citer = _cal_energy_5x5.find(layer);
  if (citer == _cal_energy_5x5.end()) return NAN;
  return citer->second;
}

unsigned int SvtxTrack_v1::get_cal_cluster_id(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, int>::const_iterator citer = _cal_cluster_id.find(layer);
  if (citer == _cal_cluster_id.end()) return -9999;
  return citer->second;
}

TrkrDefs::cluskey SvtxTrack_v1::get_cal_cluster_key(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, TrkrDefs::cluskey>::const_iterator citer = _cal_cluster_key.find(layer);
  if (citer == _cal_cluster_key.end()) return -9999;
  return citer->second;
}

float SvtxTrack_v1::get_cal_cluster_e(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, float>::const_iterator citer = _cal_cluster_e.find(layer);
  if (citer == _cal_cluster_e.end()) return NAN;
  return citer->second;
}
