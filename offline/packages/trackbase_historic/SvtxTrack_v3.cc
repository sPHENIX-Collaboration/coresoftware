#include "SvtxTrack_v3.h"
#include "SvtxTrackState.h"
#include "SvtxTrackState_v1.h"

#include <trackbase/TrkrDefs.h>  // for cluskey

#include <phool/PHObject.h>      // for PHObject

#include <climits>
#include <map>
#include <vector>                // for vector

SvtxTrack_v3::SvtxTrack_v3()
{
  for(int i = 0; i < 6; i++)
    for(int j = 0; j < 6; j++)
  { _acts_trajectory_covariance[i][j] = NAN; }

  // always include the pca point
  _states.insert( std::make_pair(0, new SvtxTrackState_v1(0)));

}

SvtxTrack_v3::SvtxTrack_v3(const SvtxTrack& source)
{ SvtxTrack_v3::CopyFrom( source ); }

// have to suppress uninitMenberVar from cppcheck since it triggers many false positive
// cppcheck-suppress uninitMemberVar
SvtxTrack_v3::SvtxTrack_v3(const SvtxTrack_v3& source)
{ SvtxTrack_v3::CopyFrom( source ); }

SvtxTrack_v3& SvtxTrack_v3::operator=(const SvtxTrack_v3& source)
{ if( this != &source ) CopyFrom( source ); return *this; }

SvtxTrack_v3::~SvtxTrack_v3()
{ clear_states(); }

void SvtxTrack_v3::CopyFrom( const SvtxTrack& source )
{
  // do nothing if copying onto oneself
  if( this == &source ) return;
  
  // parent class method
  SvtxTrack::CopyFrom( source );
  
  // copy acts covariance 
  for( int i = 0; i<6; ++i )
    for( int j = 0; j<6; ++j )
  { set_acts_covariance( i, j, source.get_acts_covariance( i, j ) ); }

  _track_id = source.get_id();
  _vertex_id = source.get_vertex_id();
  _is_positive_charge = source.get_positive_charge();
  _chisq = source.get_chisq();
  _ndf = source.get_ndf();
  _track_crossing = source.get_crossing();
  _dca = source.get_dca();
  _dca_error = source.get_dca_error();
  _dca2d = source.get_dca2d();
  _dca2d_error = source.get_dca2d_error();
  _dca3d_xy = source.get_dca3d_xy();
  _dca3d_xy_error = source.get_dca3d_xy_error();
  _dca3d_z = source.get_dca3d_z();
  _dca3d_z_error = source.get_dca3d_z_error();
  
  // copy the states over into new state objects stored here
  clear_states();
  for( auto iter = source.begin_states(); iter != source.end_states(); ++iter )
  { _states.insert( std::make_pair(iter->first, static_cast<SvtxTrackState*>(iter->second->CloneMe() ) ) ); }

  // copy over cluster ID set
  _cluster_ids.clear();
  std::copy( source.begin_clusters(), source.end_clusters(), std::inserter( _cluster_ids, _cluster_ids.begin() ) );
  
  // copy over cluster key set
  _cluster_keys.clear();
  std::copy( source.begin_cluster_keys(), source.end_cluster_keys(), std::inserter( _cluster_keys, _cluster_keys.begin() ) );

  // copy over calorimeter projections
  _cal_dphi.clear();
  _cal_deta.clear();
  _cal_energy_3x3.clear();
  _cal_energy_5x5.clear();
  _cal_cluster_id.clear();
  _cal_cluster_key.clear();
  _cal_cluster_e.clear();

  for( const auto& type: { SvtxTrack::PRES, SvtxTrack::CEMC, SvtxTrack::HCALIN, SvtxTrack::HCALOUT } )
  {
    if(!std::isnan(source.get_cal_dphi(type))) set_cal_dphi(type, source.get_cal_dphi(type));
    if(!std::isnan(source.get_cal_deta(type))) set_cal_deta(type, source.get_cal_deta(type));
    if(!std::isnan(source.get_cal_energy_3x3(type))) set_cal_energy_3x3(type, source.get_cal_energy_3x3(type));
    if(!std::isnan(source.get_cal_energy_5x5(type))) set_cal_energy_5x5(type, source.get_cal_energy_5x5(type));
    if(source.get_cal_cluster_id(type) != UINT_MAX) set_cal_cluster_id(type, source.get_cal_cluster_id(type));
    if(source.get_cal_cluster_key(type) != UINT_MAX) set_cal_cluster_key(type, source.get_cal_cluster_key(type));
    if(!std::isnan(source.get_cal_cluster_e(type))) set_cal_cluster_e(type, source.get_cal_cluster_e(type));
  }

}

void SvtxTrack_v3::identify(std::ostream& os) const
{
  os << "SvtxTrack_v3 Object ";
  os << "id: " << get_id() << " ";
  os << "vertex id: " << get_vertex_id() << " ";
  os << "charge: " << get_charge() << " ";
  os << "chisq: " << get_chisq() << " ndf:" << get_ndf() << " ";
  os << std::endl;

  os << "(px,py,pz) = ("
     << get_px() << ","
     << get_py() << ","
     << get_pz() << ")" << std::endl;

  os << "(x,y,z) = (" << get_x() << "," << get_y() << "," << get_z() << ")" << std::endl;

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
    os << " track has no clusters " << std::endl;
  
  os << std::endl;

  return;
}

void SvtxTrack_v3::clear_states()
{
  for( const auto& pair:_states )
  { delete pair.second; }
  
  _states.clear();
}

int SvtxTrack_v3::isValid() const
{
  return 1;
}

const SvtxTrackState* SvtxTrack_v3::get_state(float pathlength) const
{
  const auto iter = _states.find(pathlength);
  return (iter == _states.end()) ? nullptr:iter->second;
}

SvtxTrackState* SvtxTrack_v3::get_state(float pathlength)
{
  const auto iter = _states.find(pathlength);
  return (iter == _states.end()) ? nullptr:iter->second;
}

SvtxTrackState* SvtxTrack_v3::insert_state(const SvtxTrackState* state)
{
  // find closest iterator
  const auto pathlength = state->get_pathlength();
  auto iterator = _states.lower_bound( pathlength );
  if( iterator == _states.end() || pathlength < iterator->first )
  {
    // pathlength not found. Make a copy and insert
    const auto copy =  static_cast<SvtxTrackState*> (state->CloneMe());
    iterator = _states.insert(iterator, std::make_pair( pathlength, copy ));
  }

  // return matching state
  return iterator->second;
}

size_t SvtxTrack_v3::erase_state(float pathlength)
{
  StateIter iter = _states.find(pathlength);
  if (iter == _states.end()) return _states.size();

  delete iter->second;
  _states.erase(iter);
  return _states.size();
}

float SvtxTrack_v3::get_cal_dphi(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, float>::const_iterator citer = _cal_dphi.find(layer);
  if (citer == _cal_dphi.end()) return NAN;
  return citer->second;
}

float SvtxTrack_v3::get_cal_deta(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, float>::const_iterator citer = _cal_deta.find(layer);
  if (citer == _cal_deta.end()) return NAN;
  return citer->second;
}

float SvtxTrack_v3::get_cal_energy_3x3(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, float>::const_iterator citer = _cal_energy_3x3.find(layer);
  if (citer == _cal_energy_3x3.end()) return NAN;
  return citer->second;
}

float SvtxTrack_v3::get_cal_energy_5x5(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, float>::const_iterator citer = _cal_energy_5x5.find(layer);
  if (citer == _cal_energy_5x5.end()) return NAN;
  return citer->second;
}

unsigned int SvtxTrack_v3::get_cal_cluster_id(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, int>::const_iterator citer = _cal_cluster_id.find(layer);
  if (citer == _cal_cluster_id.end()) return UINT_MAX;
  return citer->second;
}

TrkrDefs::cluskey SvtxTrack_v3::get_cal_cluster_key(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, TrkrDefs::cluskey>::const_iterator citer = _cal_cluster_key.find(layer);
  if (citer == _cal_cluster_key.end()) return UINT_MAX;
  return citer->second;
}

float SvtxTrack_v3::get_cal_cluster_e(SvtxTrack::CAL_LAYER layer) const
{
  std::map<SvtxTrack::CAL_LAYER, float>::const_iterator citer = _cal_cluster_e.find(layer);
  if (citer == _cal_cluster_e.end()) return NAN;
  return citer->second;
}


