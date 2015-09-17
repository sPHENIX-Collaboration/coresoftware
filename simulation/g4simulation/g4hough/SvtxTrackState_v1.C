#include "SvtxTrackState_v1.h"

ClassImp(SvtxTrackState_v1)

using namespace std;

SvtxTrackState_v1::SvtxTracState_v1(float pathlength)
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

float SvtxTrackState::get_error(unsigned int i, unsigned int j) const {
  return _covar[covar_index(i,j)];
}

void SvtxTrackState::set_error(unsigned int i, unsigned int j, float value) {
  _covar[covar_index(i,j)] = value;
  return;
}

unsigned int SvtxTrackState::covar_index(unsigned int i, unsigned int j) const {
  if (i>j) std::swap(i,j);
  return i+1+(j+1)*(j)/2-1;
}
