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
  _states.insert(make_pair(0.0,State(0.0)));
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

int SvtxTrack_v1::isValid() const {
  return 1;
}

const SvtxTrack_v1::State* SvtxTrack_v1::get_state(float pathlength) const {
  ConstStateIter iter = _states.find(pathlength);
  if (iter == _states.end()) return NULL;  
  return &iter->second;
}

SvtxTrack_v1::State* SvtxTrack_v1::get_state(float pathlength) {
  StateIter iter = _states.find(pathlength);
  if (iter == _states.end()) return NULL;
  return &iter->second;
}

SvtxTrack_v1::State* SvtxTrack_v1::insert_state(const State &state) {
  float pathlength = state.get_pathlength();
  _states.insert(make_pair( pathlength , SvtxTrack_v1::State(state) ));
  return (&_states[pathlength]);
}

float SvtxTrack_v1::get_cal_energy_3x3(SvtxTrack::CAL_LAYER layer) const {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) return NAN;
  return citer->second.get_energy_3x3();
}

void SvtxTrack_v1::set_cal_energy_3x3(SvtxTrack::CAL_LAYER layer, float energy_3x3) {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) {
    _calo_matches[layer] = SvtxTrack::CaloProjection_v1();
  }
  _calo_matches[layer].set_energy_3x3(energy_3x3);
  return;
}

unsigned int SvtxTrack_v1::get_cal_cluster_id(SvtxTrack::CAL_LAYER layer) const {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) return UINT_MAX;
  return citer->second.get_cluster_id();
}

void SvtxTrack_v1::set_cal_cluster_id(SvtxTrack::CAL_LAYER layer, unsigned int clus_id) {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) {
    _calo_matches[layer] = SvtxTrack::CaloProjection_v1();
  }
  _calo_matches[layer].set_cluster_id(clus_id);
  return;
}

float SvtxTrack_v1::get_cal_dphi(SvtxTrack::CAL_LAYER layer) const {  
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) return NAN;
  return citer->second.get_dphi();
}

void SvtxTrack_v1::set_cal_dphi(SvtxTrack::CAL_LAYER layer, float dphi) {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) {
    _calo_matches[layer] = SvtxTrack::CaloProjection_v1();
  }
  _calo_matches[layer].set_dphi(dphi);
  return;
}

float SvtxTrack_v1::get_cal_deta(SvtxTrack::CAL_LAYER layer) const {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) return NAN;
  return citer->second.get_deta();
}

void SvtxTrack_v1::set_cal_deta(SvtxTrack::CAL_LAYER layer, float deta) {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) {
    _calo_matches[layer] = SvtxTrack::CaloProjection();
  }
  _calo_matches[layer].set_deta(deta);
  return;
}

float SvtxTrack_v1::get_cal_cluster_e(SvtxTrack::CAL_LAYER layer) const {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) return NAN;
  return citer->second.get_cluster_energy();
}

void SvtxTrack_v1::set_cal_cluster_e(SvtxTrack::CAL_LAYER layer, float clus_e) {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) {
    _calo_matches[layer] = SvtxTrack::CaloProjection_v1();
  }
  _calo_matches[layer].set_cluster_energy(clus_e);
  return;
}

// --- innner State class ----------------------------------------------------//

SvtxTrack_v1::State_v1::State_v1(float pathlength)
  : _pathlength(pathlength),
    _pos(),
    _mom(),
    _covar() {
  for (int i=0;i<3;++i) _pos[i] = 0.0;
  for (int i=0;i<3;++i) _mom[i] = NAN;
  for (int i = 0; i < 6; ++i) {
    for (int j = i; j < 6; ++j) {
      set_error(i,j,0.0);
    }
  } 
}

float SvtxTrack_v1::State_v1::get_error(unsigned int i, unsigned int j) const {
  return _covar[covar_index(i,j)];
}

void SvtxTrack_v1::State_v1::set_error(unsigned int i, unsigned int j, float value) {
  _covar[covar_index(i,j)] = value;
  return;
}

unsigned int SvtxTrack_v1::State_v1::covar_index(unsigned int i, unsigned int j) const {
  if (i>j) std::swap(i,j);
  return i+1+(j+1)*(j)/2-1;
}

// --- innner CaloProjection class -------------------------------------------//

SvtxTrack_v1::CaloProjection_v1::CaloProjection_v1()
  : _e3x3(NAN),
    _clus_id(UINT_MAX),
    _deta(NAN),
    _dphi(NAN),
    _clus_e(NAN) {
}
