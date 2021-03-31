#include "SvtxTrack_v2.h"
#include "SvtxTrackState.h"
#include "SvtxTrackState_v1.h"

#include <trackbase/TrkrDefs.h>  // for cluskey

#include <phool/PHObject.h>      // for PHObject

#include <climits>
#include <map>
#include <vector>                // for vector


using namespace std;

SvtxTrack_v2::SvtxTrack_v2()
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
  , _acts_mj(nullptr)
{
  // always include the pca point
  _states.insert(make_pair(0.0, new SvtxTrackState_v1(0.0)));
}

SvtxTrack_v2::SvtxTrack_v2(const SvtxTrack_v2& track)
{
  *this = track;
  return;
}

SvtxTrack_v2& SvtxTrack_v2::operator=(const SvtxTrack_v2& track)
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

SvtxTrack_v2::~SvtxTrack_v2()
{
  clear_states();
}

void SvtxTrack_v2::identify(std::ostream& os) const
{
  os << "SvtxTrack_v2 Object ";
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

void SvtxTrack_v2::clear_states()
{
  while(_states.begin() != _states.end())
  {
    delete _states.begin()->second;
    _states.erase(_states.begin());
  }
}

int SvtxTrack_v2::isValid() const
{
  return 1;
}

const SvtxTrackState* SvtxTrack_v2::get_state(float pathlength) const
{
  ConstStateIter iter = _states.find(pathlength);
  if (iter == _states.end()) return nullptr;
  return iter->second;
}

SvtxTrackState* SvtxTrack_v2::get_state(float pathlength)
{
  StateIter iter = _states.find(pathlength);
  if (iter == _states.end()) return nullptr;
  return iter->second;
}

SvtxTrackState* SvtxTrack_v2::insert_state(const SvtxTrackState* state)
{
  _states.insert(make_pair(state->get_pathlength(), dynamic_cast< SvtxTrackState*> (state->CloneMe())));
  return _states[state->get_pathlength()];
}

size_t SvtxTrack_v2::erase_state(float pathlength)
{
  StateIter iter = _states.find(pathlength);
  if (iter == _states.end()) return _states.size();

  delete iter->second;
  _states.erase(iter);
  return _states.size();
}

float SvtxTrack_v2::get_cal_dphi(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, float>::const_iterator citer = _cal_dphi.find(layer);
  if (citer == _cal_dphi.end()) return NAN;
  return citer->second;
}

float SvtxTrack_v2::get_cal_deta(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, float>::const_iterator citer = _cal_deta.find(layer);
  if (citer == _cal_deta.end()) return NAN;
  return citer->second;
}

float SvtxTrack_v2::get_cal_energy_3x3(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, float>::const_iterator citer = _cal_energy_3x3.find(layer);
  if (citer == _cal_energy_3x3.end()) return NAN;
  return citer->second;
}

float SvtxTrack_v2::get_cal_energy_5x5(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, float>::const_iterator citer = _cal_energy_5x5.find(layer);
  if (citer == _cal_energy_5x5.end()) return NAN;
  return citer->second;
}

unsigned int SvtxTrack_v2::get_cal_cluster_id(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, int>::const_iterator citer = _cal_cluster_id.find(layer);
  if (citer == _cal_cluster_id.end()) return -9999;
  return citer->second;
}

TrkrDefs::cluskey SvtxTrack_v2::get_cal_cluster_key(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, TrkrDefs::cluskey>::const_iterator citer = _cal_cluster_key.find(layer);
  if (citer == _cal_cluster_key.end()) return -9999;
  return citer->second;
}

float SvtxTrack_v2::get_cal_cluster_e(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, float>::const_iterator citer = _cal_cluster_e.find(layer);
  if (citer == _cal_cluster_e.end()) return NAN;
  return citer->second;
}


ActsTrackParametersPtr SvtxTrack_v2::get_acts_track_parameters() const
{
  Acts::Vector4D position(get_x() * Acts::UnitConstants::cm,
			  get_y() * Acts::UnitConstants::cm,
			  get_z() * Acts::UnitConstants::cm,
			  10 * Acts::UnitConstants::ns);
  
  Acts::Vector3D momentum(get_px(), get_py(), get_pz());
  double mom = get_p();
  int charge = get_charge();

  const Acts::BoundSymMatrix cov = rotateSvtxTrackCovToActs();

  return std::make_shared<ActsExamples::TrackParameters>(position, momentum,
				       mom, charge, cov);

}

Acts::BoundSymMatrix SvtxTrack_v2::rotateSvtxTrackCovToActs() const
{
  Acts::BoundSymMatrix svtxCovariance = Acts::BoundSymMatrix::Zero();
  
  for(int i = 0; i < 6; i++) {
    for(int j = 0; j < 6; j++) {
      svtxCovariance(i,j) = get_error(i,j);
      
      /// Convert Svtx to mm and GeV units as Acts expects
      if(i < 3 && j < 3)
	svtxCovariance(i,j) *= Acts::UnitConstants::cm2;
      else if (i < 3)
	svtxCovariance(i,j) *= Acts::UnitConstants::cm;
      else if (j < 3)
	svtxCovariance(i,j) *= Acts::UnitConstants::cm;
    }
  }

  double p = get_p();
  double uPx = get_px() / p;
  double uPy = get_py() / p;
  double uPz = get_pz() / p;
  
  double cosTheta = uPz;
  double sinTheta = sqrt(uPx * uPx + uPy * uPy);
  double invSinTheta = 1. / sinTheta;
  double cosPhi = uPx * invSinTheta; // equivalent to x/r
  double sinPhi = uPy * invSinTheta; // equivalent to y/r
  
  /// First we rotate to (x,y,z,time,Tx,Ty,Tz,q/p) to take advantage of the
  /// already created Acts rotation matrix from this basis into the Acts local basis
  /// We basically go backwards from rotateActsCovToSvtxTrack to get the Acts cov 
  /// from the SvtxTrack cov

  /// This is going from Acts->Svtx, so we will take the transpose
  Acts::ActsMatrixD<6,8> sphenixRot;
  sphenixRot.setZero();
  /// Make the xyz transform unity
  sphenixRot(0,0) = 1;
  sphenixRot(1,1) = 1;
  sphenixRot(2,2) = 1;
  sphenixRot(3,4) = 1./p;
  sphenixRot(4,5) = 1./p;
  sphenixRot(5,6) = 1./p;
  sphenixRot(3,7) = uPx * p * p;
  sphenixRot(4,7) = uPy * p * p;
  sphenixRot(5,7) = uPz * p * p;

  auto rotatedMatrix 
    = sphenixRot.transpose() * svtxCovariance * sphenixRot;
  
  /// Now take the 8x8 matrix and rotate it to Acts basis
  Acts::BoundToFreeMatrix jacobianLocalToGlobal = Acts::BoundToFreeMatrix::Zero();
  jacobianLocalToGlobal(0, Acts::eBoundLoc0) = -sinPhi;
  jacobianLocalToGlobal(0, Acts::eBoundLoc1) = -cosPhi * cosTheta;
  jacobianLocalToGlobal(1, Acts::eBoundLoc0) = cosPhi;
  jacobianLocalToGlobal(1, Acts::eBoundLoc1) = -sinPhi * cosTheta;
  jacobianLocalToGlobal(2, Acts::eBoundLoc1) = sinTheta;
  jacobianLocalToGlobal(3, Acts::eBoundTime) = 1;
  jacobianLocalToGlobal(4, Acts::eBoundPhi) = -sinTheta * sinPhi;
  jacobianLocalToGlobal(4, Acts::eBoundTheta) = cosTheta * cosPhi;
  jacobianLocalToGlobal(5, Acts::eBoundPhi) = sinTheta * cosPhi;
  jacobianLocalToGlobal(5, Acts::eBoundTheta) = cosTheta * sinPhi;
  jacobianLocalToGlobal(6, Acts::eBoundTheta) = -sinTheta;
  jacobianLocalToGlobal(7, Acts::eBoundQOverP) = 1;

  /// Since we are using the local to global jacobian we do R^TCR instead of 
  /// RCR^T
  Acts::BoundSymMatrix actsLocalCov = 
  jacobianLocalToGlobal.transpose() * rotatedMatrix * jacobianLocalToGlobal;

  return actsLocalCov;

}
