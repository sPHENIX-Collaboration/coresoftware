#include "SvtxTrack.h"

#include <math.h>
#include <limits.h>
#include <algorithm>
#include <map>

ClassImp(SvtxTrack)

using namespace std;

SvtxTrack::SvtxTrack()
  : _track_id(UINT_MAX),
    _is_positive_charge(false),
    _chisq(NAN),
    _ndf(0),
    _dca(NAN),
    _dca2d(NAN),
    _dca2d_error(NAN),
    _states(),
    _cluster_ids(),
    _cluster_positions(),
    _calo_matches() {
  // always include the pca point
  _states.insert(make_pair(0.0,State()));
}

void SvtxTrack::identify(std::ostream& os) const {
  os << "SvtxTrack Object ";
  os << "id: " << get_id() << " ";
  os << "charge: " << get_charge() << " ";
  os << "chisq: " << get_chisq() << " ndf:" << get_ndf() << " ";
  os << endl;

  os << "(px,py,pz) = ("
     << get_px() << ","
     << get_py() << ","
     << get_pz() << ")" << endl;

  os << "(x,y,z) = (" << get_x() << "," << get_y() << "," << get_z() << ")" << endl;
  
  if (getNhits() > 0) {
    os << "clusters: ";
    for (unsigned int i = 0; i < 100; ++i) {
      if (hasCluster(i)) {
	os << getClusterID(i) << " ";
      }
    }
  }
  os << endl;    
 
  return;
}

void SvtxTrack::Reset() {
  *this = SvtxTrack();
  return;
}

int SvtxTrack::isValid() const {
  return 1;
}

float SvtxTrack::getHitPosition(int layer, int coor) const {
  std::map<int,std::vector<float> >::const_iterator citer = _cluster_positions.find(layer);
  if (citer == _cluster_positions.end()) return NAN;
  return citer->second[coor];
}

float SvtxTrack::getInnerMostHitPosition(int coor) const {
  if (_cluster_positions.empty()) return NAN;  
  return _cluster_positions.begin()->second[coor];
}

short SvtxTrack::getNhits() const {
  return _cluster_ids.size();
}

float SvtxTrack::get_cal_energy_3x3(SvtxTrack::CAL_LAYER layer) const {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) return NAN;
  return citer->second.get_energy_3x3();
}

void SvtxTrack::set_cal_energy_3x3(SvtxTrack::CAL_LAYER layer, float energy_3x3) {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) {
    _calo_matches[layer] = SvtxTrack::CaloProjection();
  }
  _calo_matches[layer].set_energy_3x3(energy_3x3);
  return;
}

unsigned int SvtxTrack::get_cal_cluster_id(SvtxTrack::CAL_LAYER layer) const {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) return UINT_MAX;
  return citer->second.get_cluster_id();
}

void SvtxTrack::set_cal_cluster_id(SvtxTrack::CAL_LAYER layer, unsigned int clus_id) {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) {
    _calo_matches[layer] = SvtxTrack::CaloProjection();
  }
  _calo_matches[layer].set_cluster_id(clus_id);
  return;
}

float SvtxTrack::get_cal_dphi(SvtxTrack::CAL_LAYER layer) const {  
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) return NAN;
  return citer->second.get_dphi();
}

void SvtxTrack::set_cal_dphi(SvtxTrack::CAL_LAYER layer, float dphi) {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) {
    _calo_matches[layer] = SvtxTrack::CaloProjection();
  }
  _calo_matches[layer].set_dphi(dphi);
  return;
}

float SvtxTrack::get_cal_deta(SvtxTrack::CAL_LAYER layer) const {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) return NAN;
  return citer->second.get_deta();
}

void SvtxTrack::set_cal_deta(SvtxTrack::CAL_LAYER layer, float deta) {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) {
    _calo_matches[layer] = SvtxTrack::CaloProjection();
  }
  _calo_matches[layer].set_deta(deta);
  return;
}

float SvtxTrack::get_cal_cluster_e(SvtxTrack::CAL_LAYER layer) const {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) return NAN;
  return citer->second.get_cluster_energy();
}

void SvtxTrack::set_cal_cluster_e(SvtxTrack::CAL_LAYER layer, float clus_e) {
  std::map<SvtxTrack::CAL_LAYER,
	   SvtxTrack::CaloProjection>::const_iterator citer = _calo_matches.find(layer);
  if (citer == _calo_matches.end()) {
    _calo_matches[layer] = SvtxTrack::CaloProjection();
  }
  _calo_matches[layer].set_cluster_energy(clus_e);
  return;
}

// --- innner State class ----------------------------------------------------//

SvtxTrack::State::State()
  : _pos(),
    _mom(),
    _covar(),
    _helix_phi(0.0),
    _helix_d(0.0),
    _helix_kappa(0.0),
    _helix_z0(0.0),
    _helix_dzdl(0.0) {
  for (int i=0;i<3;++i) _pos[i] = 0.0;
  for (int i=0;i<3;++i) _mom[i] = NAN;
  for (int i = 0; i < 6; ++i) {
    for (int j = i; j < 6; ++j) {
      set_error(i,j,0.0);
    }
  } 
}

float SvtxTrack::State::get_error(unsigned int i, unsigned int j) const {
  return _covar[covar_index(i,j)];
}

void SvtxTrack::State::set_error(unsigned int i, unsigned int j, float value) {
  _covar[covar_index(i,j)] = value;
  return;
}

unsigned int SvtxTrack::State::covar_index(unsigned int i, unsigned int j) const {
  if (i>j) std::swap(i,j);
  return i+1+(j+1)*(j)/2-1;
}

// --- innner CaloProjection class -------------------------------------------//

SvtxTrack::CaloProjection::CaloProjection()
  : _e3x3(NAN),
    _clus_id(UINT_MAX),
    _deta(NAN),
    _dphi(NAN),
    _clus_e(NAN) {
}
