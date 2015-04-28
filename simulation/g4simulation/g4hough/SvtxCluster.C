#include "SvtxCluster.h"

#include <cmath>

#include <TMatrixF.h>

using namespace std;

ClassImp(SvtxCluster);

SvtxCluster::SvtxCluster()
  : _id(0xFFFFFFFF),
    _layer(0xFFFFFFFF),
    _pos(),
    _e(NAN),
    _adc(0xFFFFFFFF),
    _size(),
    _err(),
    _hit_ids()
{
  for (int i = 0; i < 3; ++i) _pos[i] = NAN;
  for (int i = 0; i < 3; ++i) _size[i] = new float[i+1];
  for (int i = 0; i < 3; ++i) _err[i] = new float[i+1];

  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      set_size(i,j,NAN);
      set_error(i,j,NAN);
    }
  } 
}

SvtxCluster::SvtxCluster(const SvtxCluster &clus) :
  _id(clus.get_id()),
  _layer(clus.get_layer()),
  _pos(),
  _e(clus.get_e()),
  _adc(clus.get_adc()),
  _size(),
  _err(),
  _hit_ids()
{
  for (int i=0; i<3; ++i) _pos[i] = clus.get_position(i);    
  for (int i=0; i<3; ++i) _size[i] = new float[i+1];
  for (int i=0; i<3; ++i) _err[i] = new float[i+1];
  for (int j=0; j<3; ++j) {
    for (int i=j; i<3; ++i) {
      set_size(i,j,clus.get_size(i,j));
      set_error(i,j,clus.get_error(i,j));
    }
  } 

  for (ConstHitIter iter = clus.begin_hits(); iter != clus.end_hits(); ++iter) {
    insert_hit(*iter);
  }
}

SvtxCluster& SvtxCluster::operator=(const SvtxCluster& clus) {
  Reset();

  _id = clus.get_id();
  _layer = clus.get_layer();
  for (int i=0; i<3; ++i) _pos[i] = clus.get_position(i);

  _e = clus.get_e();
  _adc = clus.get_adc();
  
  for (int j=0; j<3; ++j) {
    for (int i=j; i<3; ++i) {
      set_size(i,j,clus.get_size(i,j));
      set_error(i,j,clus.get_error(i,j));
    }
  } 

  for (ConstHitIter iter = clus.begin_hits(); iter != clus.end_hits(); ++iter) {
    insert_hit(*iter);
  }
  
  return *this;
}
  
SvtxCluster::~SvtxCluster(){
  for (int i=0; i<3; ++i) delete[] _size[i];
  for (int i=0; i<3; ++i) delete[] _err[i];
}

void SvtxCluster::identify(ostream& os) const {
  os << "---SvtxCluster-----------------------" << endl;
  os << "clusid: " << get_id() << " layer: "<< get_layer() << endl;

  os << " (x,y,z) =  (" << get_position(0);
  os << ", " << get_position(1) << ", ";
  os << get_position(2) << ") cm" << endl;

  os << " e = " << get_e() << " adc = " << get_adc() << endl;
  
  os << " size phi = " << get_phi_size();
  os << " cm, size z = " << get_z_size() << " cm" << endl;

  os << "         ( ";
  os << get_size(0,0) << " , ";
  os << get_size(0,1) << " , ";
  os << get_size(0,2) << " )" << endl;
  os << "  size = ( ";
  os << get_size(1,0) << " , ";
  os << get_size(1,1) << " , ";
  os << get_size(1,2) << " )" << endl;
  os << "         ( ";
  os << get_size(2,0) << " , ";
  os << get_size(2,1) << " , ";
  os << get_size(2,2) << " )" << endl;

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

  os << " list of hits ids: ";
  for (ConstHitIter iter = begin_hits(); iter != end_hits(); ++iter) {
    os << *iter << " ";
  }
  os << endl; 
  os << "-----------------------------------------------" << endl;
  
  return;  
}

void SvtxCluster::Reset() {
  _id    = 0xFFFFFFFF;
  _layer = 0xFFFFFFFF;
  for (int i = 0; i < 3; ++i) _pos[i] = NAN;
  _e = NAN;
  _adc = 0xFFFFFFFF;
  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      set_size(i,j,NAN);
      set_error(i,j,NAN);
    }
  } 
  _hit_ids.clear();  
}

int SvtxCluster::IsValid() const {
  if (_id == 0xFFFFFFFF) return 0;
  if (_layer == 0xFFFFFFFF) return 0;
  for (int i = 0; i < 3; ++i) {
    if (isnan(_pos[i])) return 0;
  }
  if (isnan(_e)) return 0;
  if (_adc == 0xFFFFFFFF) return 0;
  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      if (isnan(get_size(i,j))) return 0;
      if (isnan(get_error(i,j))) return 0;
    }
  }
  if (_hit_ids.empty()) return 0;

  return 1;
}

void SvtxCluster::set_size(int i, int j, float value) {
  if (j > i) set_size(j,i,value);
  else _size[i][j] = value;
  return;
}

float SvtxCluster::get_size(int i, int j) const {
  if (j > i) return get_size(j,i);
  return _size[i][j];
}

void SvtxCluster::set_error(int i, int j, float value) {
  if (j > i) set_error(j,i,value);
  else _err[i][j] = value;
  return;
}

float SvtxCluster::get_error(int i, int j) const {
  if (j > i) return get_error(j,i);
  return _err[i][j];
}

float SvtxCluster::get_phi_size() const {

  TMatrixF COVAR(3,3);
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      COVAR[i][j] = get_size(i,j);
    }
  }

  float phi = -1.0*atan2(_pos[1],_pos[0]);
  
  TMatrixF ROT(3,3);
  ROT[0][0] = cos(phi);
  ROT[0][1] = -sin(phi);
  ROT[0][2] = 0.0;
  ROT[1][0] = sin(phi);
  ROT[1][1] = cos(phi);
  ROT[1][2] = 0.0;
  ROT[2][0] = 0.0;
  ROT[2][1] = 0.0;
  ROT[2][2] = 1.0;

  TMatrixF ROT_T(3,3);
  ROT_T.Transpose(ROT);
  
  TMatrixF TRANS(3,3);
  TRANS = ROT * COVAR * ROT_T;
  
  return 2.0*sqrt(TRANS[1][1]);
}

float SvtxCluster::get_z_size() const {
  return 2.0*sqrt(get_size(2,2));
}
