#include "SvtxVertex.h"

#include <cmath>

using namespace std;

ClassImp(SvtxVertex);

SvtxVertex::SvtxVertex()
  : _id(0xFFFFFFFF),
    _t0(NAN),
    _pos(),
    _chisq(NAN),
    _ndof(0xFFFFFFFF),
    _err(3),
    _track_ids() {
  
  for (int i = 0; i < 3; ++i) _pos[i] = NAN;  
  for (int i = 0; i < 3; ++i) _err[i] = std::vector<float>(i+1);

  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      set_error(i,j,NAN);
    }
  } 
}

SvtxVertex::SvtxVertex(const SvtxVertex &vertex) :
  _id(vertex.get_id()),
  _t0(vertex.get_t0()),
  _pos(),
  _chisq(vertex.get_chisq()),
  _ndof(vertex.get_ndof()),
  _err(3),
  _track_ids() {
  
  for (int i=0; i<3; ++i) _pos[i] = vertex.get_position(i);    
  for (int i=0; i<3; ++i) _err[i] = std::vector<float>(i+1);
  for (int j=0; j<3; ++j) {
    for (int i=j; i<3; ++i) {
      set_error(i,j,vertex.get_error(i,j));
    }
  } 

  for (ConstTrackIter iter = vertex.begin_tracks(); iter != vertex.end_tracks(); ++iter) {
    insert_track(*iter);
  }
}

SvtxVertex& SvtxVertex::operator=(const SvtxVertex &vertex) {
  Reset();
  
  _id = vertex.get_id();
  _t0 = vertex.get_t0();
  for (int i=0; i<3; ++i) _pos[i] = vertex.get_position(i);    
  _chisq = vertex.get_chisq();
  _ndof = vertex.get_ndof();
  for (int j=0; j<3; ++j) {
    for (int i=j; i<3; ++i) {
      set_error(i,j,vertex.get_error(i,j));
    }
  } 

  for (ConstTrackIter iter = vertex.begin_tracks(); iter != vertex.end_tracks(); ++iter) {
    insert_track(*iter);
  }

  return *this;
}

SvtxVertex::~SvtxVertex(){}

void SvtxVertex::identify(ostream& os) const {
  os << "---SvtxVertex-----------------------" << endl;
  os << "vertexid: " << get_id() << endl;

  os << " t0 = " << get_t0() << endl;
  
  os << " (x,y,z) =  (" << get_position(0);
  os << ", " << get_position(1) << ", ";
  os << get_position(2) << ") cm" << endl;

  os << " chisq = " << get_chisq() << ", ";
  os << " ndof = " << get_ndof() << endl;

  os << "         ( ";
  os << get_error(0,0) << " , ";
  os << get_error(0,1) << " , ";
  os << get_error(0,2) << " )" << endl; 
  os << "  err  = ( ";
  os << get_error(1,0) << " , ";
  os << get_error(1,1) << " , ";
  os << get_error(1,2) << " )" << endl;
  os << "         ( ";
  os << get_error(2,0) << " , ";
  os << get_error(2,1) << " , ";
  os << get_error(2,2) << " )" << endl;

  os << " list of tracks ids: ";
  for (ConstTrackIter iter = begin_tracks(); iter != end_tracks(); ++iter) {
    os << *iter << " ";
  }
  os << endl; 
  os << "-----------------------------------------------" << endl;
  
  return;  
}

void SvtxVertex::Reset() {
  _id    = 0xFFFFFFFF;
  _t0 = NAN;
  _chisq = NAN;
  _ndof = 0xFFFFFFFF;
  
  for (int i = 0; i < 3; ++i) _pos[i] = NAN;
  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      set_error(i,j,NAN);
    }
  } 
  _track_ids.clear();
}

int SvtxVertex::IsValid() const {
  if (_id == 0xFFFFFFFF) return 0;
  if (isnan(_t0)) return 0;
  if (isnan(_chisq)) return 0;
  if (_ndof == 0xFFFFFFFF) return 0;
  
  for (int i = 0; i < 3; ++i) {
    if (isnan(_pos[i])) return 0;
  }
  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      if (isnan(get_error(i,j))) return 0;
    }
  }
  if (_track_ids.empty()) return 0;
  return 1;
}

void SvtxVertex::set_error(int i, int j, float value) {
  if (j > i) set_error(j,i,value);
  else _err[i][j] = value;
  return;
}

float SvtxVertex::get_error(int i, int j) const {
  if (j > i) return get_error(j,i);
  return _err[i][j];
}
