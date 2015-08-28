#include "GlobalVertex_v1.h"

#include <cmath>

using namespace std;

ClassImp(GlobalVertex_v1);

GlobalVertex_v1::GlobalVertex_v1()
  : _id(0xFFFFFFFF),
    _t(NAN),
    _t_err(NAN),
    _pos(),
    _chisq(NAN),
    _ndof(0xFFFFFFFF),
    _err(3),
    _vtx_ids() {
  
  for (int i = 0; i < 3; ++i) _pos[i] = NAN;  
  for (int i = 0; i < 3; ++i) _err[i] = std::vector<float>(i+1);

  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      set_error(i,j,NAN);
    }
  } 
}

GlobalVertex_v1::~GlobalVertex_v1() {}

void GlobalVertex_v1::identify(ostream& os) const {
  os << "---GlobalVertex_v1-----------------------" << endl;
  os << "vertexid: " << get_id() << endl;

  os << " t = " << get_t() << endl;
  
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

  os << " list of vtx ids: ";
  for (ConstVtxIter iter = begin_vtxids(); iter != end_vtxids(); ++iter) {
    os << iter->first << " => " << iter->second << endl;
  }
  os << "-----------------------------------------------" << endl;
  
  return;  
}

void GlobalVertex_v1::Reset() {
  _id    = 0xFFFFFFFF;
  _t = NAN;
  _t_err = NAN;
  _chisq = NAN;
  _ndof = 0xFFFFFFFF;
  
  for (int i = 0; i < 3; ++i) _pos[i] = NAN;
  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      set_error(i,j,NAN);
    }
  } 
  _vtx_ids.clear();
}

int GlobalVertex_v1::IsValid() const {
  if (_id == 0xFFFFFFFF) return 0;
  if (isnan(_t)) return 0;
  if (isnan(_t_err)) return 0;
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
  if (_vtx_ids.empty()) return 0;
  return 1;
}

void GlobalVertex_v1::set_error(int i, int j, float value) {
  if (j > i) set_error(j,i,value);
  else _err[i][j] = value;
  return;
}

float GlobalVertex_v1::get_error(int i, int j) const {
  if (j > i) return get_error(j,i);
  return _err[i][j];
}
