#include "SvtxTrack_v4.h"
#include "SvtxTrackState.h"
#include "SvtxTrackState_v1.h"

#include <trackbase/TrkrDefs.h>  // for cluskey

#include <phool/PHObject.h>      // for PHObject

#include <climits>
#include <map>
#include <vector>                // for vector

SvtxTrack_v4::SvtxTrack_v4()
{
  // always include the pca point
  _states.insert( std::make_pair(0, new SvtxTrackState_v1(0)));

}

SvtxTrack_v4::SvtxTrack_v4(const SvtxTrack& source)
{ SvtxTrack_v4::CopyFrom( source ); }

// have to suppress uninitMenberVar from cppcheck since it triggers many false positive
// cppcheck-suppress uninitMemberVar
SvtxTrack_v4::SvtxTrack_v4(const SvtxTrack_v4& source)
{ SvtxTrack_v4::CopyFrom( source ); }

SvtxTrack_v4& SvtxTrack_v4::operator=(const SvtxTrack_v4& source)
{ CopyFrom( source ); return *this; }

SvtxTrack_v4::~SvtxTrack_v4()
{ clear_states(); }

void SvtxTrack_v4::CopyFrom( const SvtxTrack& source )
{
  // do nothing if copying onto oneself
  if( this == &source ) return;
  
  // parent class method
  SvtxTrack::CopyFrom( source );
  
  _tpc_seed = source.get_tpc_seed();
  _silicon_seed = source.get_silicon_seed();
  _vertex_id = source.get_vertex_id();
  _is_positive_charge = source.get_positive_charge();
  _chisq = source.get_chisq();
  _ndf = source.get_ndf();
  _track_crossing = source.get_crossing();
  
  // copy the states over into new state objects stored here
  clear_states();
  for( auto iter = source.begin_states(); iter != source.end_states(); ++iter )
  { _states.insert( std::make_pair(iter->first, static_cast<SvtxTrackState*>(iter->second->CloneMe() ) ) ); }

}

void SvtxTrack_v4::identify(std::ostream& os) const
{
  os << "SvtxTrack_v4 Object ";
  os << "id: " << get_id() << " ";
  os << "vertex id: " << get_vertex_id() << " ";
  os << "charge: " << get_charge() << " ";
  os << "chisq: " << get_chisq() << " ndf:" << get_ndf() << " ";
  os << "nstates: " << _states.size() << " ";
  os << std::endl;

  os << "(px,py,pz) = ("
     << get_px() << ","
     << get_py() << ","
     << get_pz() << ")" << std::endl;

  os << "(x,y,z) = (" << get_x() << "," << get_y() << "," << get_z() << ")" << std::endl;

  os << "Silicon clusters " << std::endl;
  if(_silicon_seed)
    {
      for(auto iter = _silicon_seed->begin_cluster_keys(); 
	  iter != _silicon_seed->end_cluster_keys();
	  ++iter)
	{
	  std::cout << *iter << ", ";
	}
    }
  os << std::endl << "Tpc + TPOT clusters " << std::endl;
  if(_tpc_seed)
  {
    for(auto iter = _tpc_seed->begin_cluster_keys(); 
      iter != _tpc_seed->end_cluster_keys();
      ++iter)
    {
      std::cout << *iter << ", ";
    }
  }
  os << std::endl;

  return;
}

void SvtxTrack_v4::clear_states()
{
  for( const auto& pair:_states )
  { delete pair.second; }
  
  _states.clear();
}

int SvtxTrack_v4::isValid() const
{
  return 1;
}

const SvtxTrackState* SvtxTrack_v4::get_state(float pathlength) const
{
  const auto iter = _states.find(pathlength);
  return (iter == _states.end()) ? nullptr:iter->second;
}

SvtxTrackState* SvtxTrack_v4::get_state(float pathlength)
{
  const auto iter = _states.find(pathlength);
  return (iter == _states.end()) ? nullptr:iter->second;
}

SvtxTrackState* SvtxTrack_v4::insert_state(const SvtxTrackState* state)
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

size_t SvtxTrack_v4::erase_state(float pathlength)
{
  StateIter iter = _states.find(pathlength);
  if (iter == _states.end()) return _states.size();

  delete iter->second;
  _states.erase(iter);
  return _states.size();
}


