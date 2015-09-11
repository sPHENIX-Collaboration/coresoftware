#ifndef __SVTXTRACK_H__
#define __SVTXTRACK_H__

#include <phool/PHObject.h>

#include <TMatrix.h>

#include <iostream>
#include <map>

class SvtxTrack : public PHObject {
  
 public:

  enum CAL_LAYER {PRES,CEMC,HCALIN,HCALOUT};

  SvtxTrack();
  virtual ~SvtxTrack() {}
  
  // The "standard PHObject response" functions...
  void identify(std::ostream &os=std::cout) const;
  void Reset();
  int  isValid() const;

  void set_id(int id) {_track_id = id;}
  void setTrackID(int index){_track_id = index;}
  int getTrackID() const {return _track_id;}  
  
  void setClusterID(int layer, int index) {_clusterID[layer] = index;}
  int getClusterID(int layer) const {return _clusterID[layer];}
  bool hasCluster(int layer) const {return (_clusterID[layer] >- 9999);}
  
  void setScatter(int layer, float sct) {_scatter[layer] = sct;}
  float getScatter(int layer) const {return _scatter[layer];}
  
  void setHitPosition(int layer, float x, float y, float z) {
    _position[layer][0]=x;
    _position[layer][1]=y;
    _position[layer][2]=z;
  }
  float getHitPosition(int layer, int coor) const {return _position[layer][coor];}

  void setMomentum(float p) {_momentum = p;}
  float getMomentum() const {return _momentum;}
  
  void set3Momentum(float px, float py, float pz) {
    _mom3[0] = px;
    _mom3[1] = py;
    _mom3[2] = pz;
  };
  float get3Momentum(int coor) const {return _mom3[coor];}
  
  void setCharge(int c) {_charge = c;}
  int getCharge() const {return _charge;}
  
  void setPrimary(bool prim) {}
  bool getPrimary() const {return false;}
  
  void setPositive(bool prim) {_ispositive = prim;}
  bool getPositive() const {return _ispositive;}
  
  //void setNhits(int layer, short n);
  short getNhits() const;
  
  void setQuality(float q) {_quality = q;}
  float getQuality() const {return _quality;}

  void setChisq(float q) {_chisq = q;}
  float getChisq() const {return _chisq;}

  void setChisqv(float q) {_chisqv = q;}
  float getChisqv() const {return _chisqv;}

  void setNDF(int q) {_ndf = q;}
  int getNDF() const {return _ndf;}

  void setDCA(float d) {_DCA = d;}
  float getDCA() const {return _DCA;}

  void setDCA2D(float d) {_DCA2D = d;}
  float getDCA2D() const {return _DCA2D;}
  
  void setDCA2Dsigma(float s) {_DCA2Dsigma = s;}
  float getDCA2Dsigma() const {return _DCA2Dsigma;}

  float getInnerMostHitPosition(int coor) const;

  void  set_phi(float phi) {_phi = phi;}
  float get_phi() const {return _phi;}

  void  set_d(float d) {_d = d;}
  float get_d() const {return _d;}
  
  void  set_kappa(float kappa) {_kappa = kappa;}
  float get_kappa() const {return _kappa;}
    
  void set_z0(float z0) {_z0 = z0;}
  float get_z0() const {return _z0;}

  void  set_dzdl(float dzdl) {_dzdl = dzdl;}
  float get_dzdl() const {return _dzdl;}
  
  const TMatrix* getCovariance() const {return &_covariance;}
  TMatrix* getCovariance() {return &_covariance;}  

  void set_cal_dphi(int layer, float dphi) {_cal_dphi[layer] = dphi;}
  float get_cal_dphi(int layer) const {return _cal_dphi[layer];}

  void set_cal_deta(int layer, float deta) {_cal_deta[layer] = deta;}
  float get_cal_deta(int layer) const {return _cal_deta[layer];}

  void set_cal_energy_3x3(int layer, float energy_3x3) {_cal_energy_3x3[layer] = energy_3x3;}
  float get_cal_energy_3x3(int layer) const {return _cal_energy_3x3[layer];}

  void set_cal_cluster_id(int layer, int id) {_cal_cluster_id[layer] = id;}
  float get_cal_cluster_id(int layer) const {return _cal_cluster_id[layer];}

  void set_cal_cluster_e(int layer, float e) {_cal_cluster_e[layer] = e;}
  float get_cal_cluster_e(int layer) const {return _cal_cluster_e[layer];}

  float get_x() const{return _x;}
  void set_x(float val){_x = val;}
  float get_y() const{return _y;}
  void set_y(float val){_y = val;}
  float get_z() const{return _z;}
  void set_z(float val){_z = val;}

 private: 

  int     _track_id;
  float   _phi,_d,_kappa,_z0,_dzdl;
  int     _clusterID[100];
  float   _position[100][3];
  float   _momentum;
  float   _mom3[3];
  int     _charge;
  bool    _ispositive;
  float   _quality;
  float   _chisq;
  float   _chisqv;
  int     _ndf;
  float   _DCA;
  float   _DCA2D;
  float   _DCA2Dsigma;
  float   _scatter[100];
  float   _x,_y,_z;
  
  TMatrix _covariance;
  
  // calorimeter matches
  float   _cal_dphi[4];
  float   _cal_deta[4];
  float   _cal_energy_3x3[4];
  int     _cal_cluster_id[4];
  float   _cal_cluster_e[4];

  // cluster ids
  //std::multimap<unsigned int,unsigned int> _cluster_ids; //< layer index => cluster id
  
  ClassDef(SvtxTrack,1)
};

#endif

