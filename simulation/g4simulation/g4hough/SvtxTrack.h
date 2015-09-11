#ifndef __SVTXTRACK_H__
#define __SVTXTRACK_H__

#include <phool/PHObject.h>

#include <TMatrix.h>

#include <iostream>
#include <map>

class SvtxTrack : public PHObject {
  
 public:

  enum CAL_LAYER {PRES=0,CEMC=1,HCALIN=2,HCALOUT=3};

  SvtxTrack();
  virtual ~SvtxTrack() {}
  
  // The "standard PHObject response" functions...
  void identify(std::ostream &os=std::cout) const;
  void Reset();
  int  isValid() const;

  void set_id(int id) {_track_id = id;}
  void setTrackID(int index){_track_id = index;}
  int getTrackID() const {return _track_id;}  
  
  void setClusterID(int layer, int index) {_cluster_ids[layer] = index;}
  int getClusterID(int layer) const {return _cluster_ids.find(layer)->second;}
  bool hasCluster(int layer) const {return (_cluster_ids.find(layer) != _cluster_ids.end());}
  
  void setScatter(int layer, float sct) {
    std::cout << "SvtxTrack:: ERROR - deprecated interface call" << std::endl;
  }
  float getScatter(int layer) const {
    std::cout << "SvtxTrack:: ERROR - deprecated interface call" << std::endl;
    return NAN;
  }
  
  void setHitPosition(int layer, float x, float y, float z) {
    std::vector<float> position(3);
    position[0] = x;
    position[1] = y;
    position[2] = z;
    _cluster_positions[layer] = position;
  }
  float getHitPosition(int layer, int coor) const;

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
  
  void setPositive(bool prim) {_ispositive = prim;}
  bool getPositive() const {return _ispositive;}

  void setPrimary(bool prim) {}
  bool getPrimary() const {return false;}
  
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

  void  set_cal_dphi(CAL_LAYER layer, float dphi) {_cal_dphi[layer] = dphi;}
  float get_cal_dphi(CAL_LAYER layer) const;

  void  set_cal_deta(CAL_LAYER layer, float deta) {_cal_deta[layer] = deta;}
  float get_cal_deta(CAL_LAYER layer) const;

  void  set_cal_energy_3x3(CAL_LAYER layer, float energy_3x3) {_cal_energy_3x3[layer] = energy_3x3;}
  float get_cal_energy_3x3(CAL_LAYER layer) const;

  void  set_cal_cluster_id(CAL_LAYER layer, int id) {_cal_cluster_id[layer] = id;}
  int   get_cal_cluster_id(CAL_LAYER layer) const;

  void  set_cal_cluster_e(CAL_LAYER layer, float e) {_cal_cluster_e[layer] = e;}
  float get_cal_cluster_e(CAL_LAYER layer) const;

  float get_x() const{return _x;}
  void set_x(float val){_x = val;}
  float get_y() const{return _y;}
  void set_y(float val){_y = val;}
  float get_z() const{return _z;}
  void set_z(float val){_z = val;}

 private: 

  int     _track_id;

  // track information
  int     _charge;
  bool    _ispositive;
  float   _quality;
  float   _chisq;
  float   _chisqv;
  int     _ndf;

  float   _DCA;
  float   _DCA2D;
  float   _DCA2Dsigma;

  // projection information
  float   _phi,_d,_kappa,_z0,_dzdl;
  float   _momentum;
  float   _mom3[3];
  float   _x,_y,_z; 
  TMatrix _covariance;
  
  // cluster contents
  std::map<int,int> _cluster_ids;                       //< layer index => cluster id
  std::map<int,std::vector<float> > _cluster_positions; //< layer index => (x,y,z)
  
  // calorimeter matches
  std::map<CAL_LAYER,float> _cal_dphi;
  std::map<CAL_LAYER,float> _cal_deta;
  std::map<CAL_LAYER,float> _cal_energy_3x3;
  std::map<CAL_LAYER,int>   _cal_cluster_id;
  std::map<CAL_LAYER,float> _cal_cluster_e;
  
  ClassDef(SvtxTrack,1)
};

#endif

