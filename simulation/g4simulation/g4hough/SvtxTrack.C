#include "SvtxTrack.h"

#include <math.h>

ClassImp(SvtxTrack)

using namespace std;

SvtxTrack::SvtxTrack()
  : _track_id(-1),
    _is_positive_charge(false),
    _quality(NAN),
    _chisq(NAN),
    _ndf(0),
    _DCA(NAN),
    _DCA2D(NAN),
    _DCA2Dsigma(NAN),
    _phi(0.0),
    _d(0.0),
    _kappa(0.0),
    _z0(0.0),
    _dzdl(0.0),
    _momentum(NAN),
    _mom3(),
    _x(0.0),
    _y(0.0),
    _z(0.0),
    _covariance(6,6),
    _cluster_ids(),
    _cluster_positions(),
    _cal_dphi(),
    _cal_deta(),
    _cal_energy_3x3(),
    _cal_cluster_id(),
    _cal_cluster_e() {
  Reset();
}

void SvtxTrack::identify(std::ostream& os) const {
  os << "SvtxTrack Object ";
  os << "id: " << getTrackID() << " ";

  if (getNhits() > 0) {
    os << "clusters: ";
    for (unsigned int i = 0; i < 100; ++i) {
      if (hasCluster(i)) {
	os << getClusterID(i) << " ";
      }
    }
  }
    
  if (getNDF() != 0) {
    os << "chisq/dof: " << getChisq()/getNDF() << " ";
  }
  os << endl;
  return;
}

void SvtxTrack::Reset() {

  _track_id = -1;

  _cluster_ids.clear();
  _cluster_positions.clear();

  _momentum=NAN;
  for(int j=0;j<3;j++){
    _mom3[j]=NAN;
  }

  _is_positive_charge = false;
  _quality=NAN;
  
  _DCA=NAN;
  _DCA2D=NAN;

  _chisq = NAN;
  _ndf = 0;
  
  _cal_dphi.clear();
  _cal_deta.clear();
  _cal_energy_3x3.clear();
  _cal_cluster_id.clear();
  _cal_cluster_e.clear();

  return;
}

int SvtxTrack::isValid() const
{
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

float SvtxTrack::get_cal_dphi(SvtxTrack::CAL_LAYER layer) const {
  std::map<SvtxTrack::CAL_LAYER,float>::const_iterator citer = _cal_dphi.find(layer);
  if (citer == _cal_dphi.end()) return NAN;
  return citer->second;
}

float SvtxTrack::get_cal_deta(SvtxTrack::CAL_LAYER layer) const {
  std::map<SvtxTrack::CAL_LAYER,float>::const_iterator citer = _cal_deta.find(layer);
  if (citer == _cal_deta.end()) return NAN;
  return citer->second;
}

float SvtxTrack::get_cal_energy_3x3(SvtxTrack::CAL_LAYER layer) const {
  std::map<SvtxTrack::CAL_LAYER,float>::const_iterator citer = _cal_energy_3x3.find(layer);
  if (citer == _cal_energy_3x3.end()) return NAN;
  return citer->second;
}

int SvtxTrack::get_cal_cluster_id(SvtxTrack::CAL_LAYER layer) const {
  std::map<SvtxTrack::CAL_LAYER,int>::const_iterator citer = _cal_cluster_id.find(layer);
  if (citer == _cal_cluster_id.end()) return -9999;
  return citer->second;
}

float SvtxTrack::get_cal_cluster_e(SvtxTrack::CAL_LAYER layer) const {
  std::map<SvtxTrack::CAL_LAYER,float>::const_iterator citer = _cal_cluster_e.find(layer);
  if (citer == _cal_cluster_e.end()) return NAN;
  return citer->second;
}
