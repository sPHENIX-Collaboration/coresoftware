#include "SvtxTrack_v1.h"

#include <math.h>
#include <limits.h>
#include <algorithm>
#include <map>

ClassImp(SvtxTrack_v1)

using namespace std;

SvtxTrack_v1::SvtxTrack_v1()
  : _track_id(UINT_MAX),
    _is_positive_charge(false),
    _chisq(NAN),
    _ndf(0),
    _dca(NAN),
    _dca2d(NAN),
    _dca2d_error(NAN),
    _states(),
    _cluster_ids(),
    _calo_matches() {
  // always include the pca point
  _states.insert(make_pair(0.0,new SvtxTrackState_v1(0.0)));
}

SvtxTrack_v1::SvtxTrack_v1(const SvtxTrack_v1& track) {
  *this = track;
  return;
}

SvtxTrack_v1& SvtxTrack_v1::operator=(const SvtxTrack_v1& track) {

  _track_id = track.get_id();
  _is_positive_charge = track.get_positive_charge();
  _chisq = track.get_chisq();
  _ndf = track.get_ndf();
  _dca = track.get_dca();
  _dca2d = track.get_dca2d();
  _dca2d_error = track.get_dca2d_error();

  // copy the states over into new state objects stored here
  clear_states();
  for (TrackStateIter iter = track.begin_states();
       iter != track.end_states();
       ++iter) {
    SvtxTrackState *state = iter->second;
    _states.insert(make_pair(state->get_pathlength(),new SvtxTrackState(state)));
  }  
  
  // copy over cluster set
  _cluster_ids.clear();
  for (ClusterConstIter iter = track.begin_clusters();
       iter != track.end_clusters();
       ++iter) {
    _cluster_ids.insert(*iter);
  }

  // copy over calorimeter projections
  _cal_dphi.clear();
  if (!isnan(track.get_cal_dphi(SvtxTrack::PRES))) set_cal_dphi(PRES,track.get_cal_dphi(SvtxTrack::PRES));
  if (!isnan(track.get_cal_dphi(SvtxTrack::CEMC))) set_cal_dphi(PRES,track.get_cal_dphi(SvtxTrack::CEMC));
  if (!isnan(track.get_cal_dphi(SvtxTrack::HCALIN))) set_cal_dphi(PRES,track.get_cal_dphi(SvtxTrack::HCALIN));
  if (!isnan(track.get_cal_dphi(SvtxTrack::HCALOUT))) set_cal_dphi(PRES,track.get_cal_dphi(SvtxTrack::HCALOUT));

  _cal_deta.clear();
  if (!isnan(track.get_cal_deta(SvtxTrack::PRES))) set_cal_deta(PRES,track.get_cal_deta(SvtxTrack::PRES));
  if (!isnan(track.get_cal_deta(SvtxTrack::CEMC))) set_cal_deta(PRES,track.get_cal_deta(SvtxTrack::CEMC));
  if (!isnan(track.get_cal_deta(SvtxTrack::HCALIN))) set_cal_deta(PRES,track.get_cal_deta(SvtxTrack::HCALIN));
  if (!isnan(track.get_cal_deta(SvtxTrack::HCALOUT))) set_cal_deta(PRES,track.get_cal_deta(SvtxTrack::HCALOUT));

  _cal_energy_3x3.clear();
  if (!isnan(track.get_cal_energy_3x3(SvtxTrack::PRES))) set_cal_energy_3x3(PRES,track.get_cal_energy_3x3(SvtxTrack::PRES));
  if (!isnan(track.get_cal_energy_3x3(SvtxTrack::CEMC))) set_cal_energy_3x3(PRES,track.get_cal_energy_3x3(SvtxTrack::CEMC));
  if (!isnan(track.get_cal_energy_3x3(SvtxTrack::HCALIN))) set_cal_energy_3x3(PRES,track.get_cal_energy_3x3(SvtxTrack::HCALIN));
  if (!isnan(track.get_cal_energy_3x3(SvtxTrack::HCALOUT))) set_cal_energy_3x3(PRES,track.get_cal_energy_3x3(SvtxTrack::HCALOUT));

  _cal_cluster_id.clear();
  if (!isnan(track.get_cal_cluster_id(SvtxTrack::PRES))) set_cal_cluster_id(PRES,track.get_cal_cluster_id(SvtxTrack::PRES));
  if (!isnan(track.get_cal_cluster_id(SvtxTrack::CEMC))) set_cal_cluster_id(PRES,track.get_cal_cluster_id(SvtxTrack::CEMC));
  if (!isnan(track.get_cal_cluster_id(SvtxTrack::HCALIN))) set_cal_cluster_id(PRES,track.get_cal_cluster_id(SvtxTrack::HCALIN));
  if (!isnan(track.get_cal_cluster_id(SvtxTrack::HCALOUT))) set_cal_cluster_id(PRES,track.get_cal_cluster_id(SvtxTrack::HCALOUT));

  _cal_cluster_e.clear();
  if (!isnan(track.get_cal_cluster_e(SvtxTrack::PRES))) set_cal_cluster_e(PRES,track.get_cal_cluster_e(SvtxTrack::PRES));
  if (!isnan(track.get_cal_cluster_e(SvtxTrack::CEMC))) set_cal_cluster_e(PRES,track.get_cal_cluster_e(SvtxTrack::CEMC));
  if (!isnan(track.get_cal_cluster_e(SvtxTrack::HCALIN))) set_cal_cluster_e(PRES,track.get_cal_cluster_e(SvtxTrack::HCALIN));
  if (!isnan(track.get_cal_cluster_e(SvtxTrack::HCALOUT))) set_cal_cluster_e(PRES,track.get_cal_cluster_e(SvtxTrack::HCALOUT));

}

SvtxTrack_v1::~SvtxTrack_v1() {
  clear_states();
}

void SvtxTrack_v1::identify(std::ostream& os) const {
  os << "SvtxTrack_v1 Object ";
  os << "id: " << get_id() << " ";
  os << "charge: " << get_charge() << " ";
  os << "chisq: " << get_chisq() << " ndf:" << get_ndf() << " ";
  os << endl;

  os << "(px,py,pz) = ("
     << get_px() << ","
     << get_py() << ","
     << get_pz() << ")" << endl;

  os << "(x,y,z) = (" << get_x() << "," << get_y() << "," << get_z() << ")" << endl;
  
  if (!empty_clusters()) {
    os << "clusters: ";
    for (SvtxTrack::ConstClusterIter iter = begin_clusters();
	 iter != end_clusters();
	 ++iter) {
      unsigned int cluster_id = *iter;
      os << cluster_id << " ";
    }
  }
  os << endl;    
 
  return;
}

void SvtxTrack_v1::clear_states() {
 for (TrackStateIter iter = _states.begin();
       iter != _state.end();
       ++iter) {
    SvtxTrackState *state = iter->second;
    delete state;
  }
  _states.clear();
}

int SvtxTrack_v1::isValid() const {
  return 1;
}

const SvtxTrack_v1::SvtxTrackState* SvtxTrack_v1::get_state(float pathlength) const {
  ConstStateIter iter = _states.find(pathlength);
  if (iter == _states.end()) return NULL;  
  return iter->second;
}

SvtxTrack_v1::SvtxTrackState* SvtxTrack_v1::get_state(float pathlength) {
  StateIter iter = _states.find(pathlength);
  if (iter == _states.end()) return NULL;
  return iter->second;
}

SvtxTrack_v1::SvtxTrackState* SvtxTrack_v1::insert_state(const SvtxTrackState* state) {
  _states.insert(make_pair( state->get_pathlength() , state->Clone() ));
  return _states[state->get_pathlength()];
}

size_t SvtxTrack_v1::erase_state(float pathlength) {
  StateIter iter = _states.find(pathlenght);
  if (iter == _states.end()) return _states.size();
  
  delete iter->second;
  return _states.erase(iter);
}
