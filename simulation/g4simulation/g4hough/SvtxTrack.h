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
  //SvtxTrack(SvtxTrack *track);
  //SvtxTrack(const SvtxTrack& track);
  virtual ~SvtxTrack() {};
  
  // The "standard PHObject response" functions...
  void identify(std::ostream &os=std::cout) const;
  void Reset();
  int  isValid() const;

  void set_id(int id) {_track_id = id;}
  void setTrackID(int index){_track_id = index;}
  int getTrackID() const {return _track_id;}  
  
  void setClusterID(int layer, int index) {clusterID[layer] = index;}
  int getClusterID(int layer) const {return clusterID[layer];}
  bool hasCluster(int layer) const {return (clusterID[layer] >- 9999);}
  
  void setScatter(int layer, float sct) {scatter[layer] = sct;}
  float getScatter(int layer) const {return scatter[layer];}
  
  void setHitPosition(int layer, float x, float y, float z) {
    position[layer][0]=x;
    position[layer][1]=y;
    position[layer][2]=z;
  }
  float getHitPosition(int layer, int coor) const {return position[layer][coor];}

  void setMomentum(float p) {momentum = p;}
  float getMomentum() const {return momentum;}
  
  void set3Momentum(float px, float py, float pz) {
    mom3[0] = px;
    mom3[1] = py;
    mom3[2] = pz;
  };
  float get3Momentum(int coor) const {return mom3[coor];}
  
  void setCharge(int c) {charge = c;}
  int getCharge() const {return charge;}
  
  void setPrimary(bool prim) {isprimary = prim;}
  bool getPrimary() const {return isprimary;}
  
  void setPositive(bool prim) {ispositive = prim;}
  bool getPositive() const {return ispositive;}
  
  //void setNhits(int layer, short n);
  short getNhits() const;
  
  void setQuality(float q) {quality = q;}
  float getQuality() const {return quality;}

  void setChisq(float q) {chisq = q;}
  float getChisq() const {return chisq;}

  void setChisqv(float q) {chisqv = q;}
  float getChisqv() const {return chisqv;}

  void setNDF(int q) {ndf = q;}
  int getNDF() const {return ndf;}

  void setDCA(float d) {DCA = d;}
  float getDCA() const {return DCA;}

  void setDCA2D(float d) {DCA2D = d;}
  float getDCA2D() const {return DCA2D;}
  
  void setDCA2Dsigma(float s) {DCA2Dsigma = s;}
  float getDCA2Dsigma() const {return DCA2Dsigma;}

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
  
  const TMatrix* getCovariance() const {return &covariance;}
  TMatrix* getCovariance() {return &covariance;}
  

  void set_cal_dphi(int layer, float dphi) {cal_dphi[layer] = dphi;}
  float get_cal_dphi(int layer) const {return cal_dphi[layer];}

  void set_cal_deta(int layer, float deta) {cal_deta[layer] = deta;}
  float get_cal_deta(int layer) const {return cal_deta[layer];}

  void set_cal_energy_3x3(int layer, float energy_3x3) {cal_energy_3x3[layer] = energy_3x3;}
  float get_cal_energy_3x3(int layer) const {return cal_energy_3x3[layer];}

  void set_cal_cluster_id(int layer, int id) {cal_cluster_id[layer] = id;}
  float get_cal_cluster_id(int layer) const {return cal_cluster_id[layer];}

  void set_cal_cluster_e(int layer, float e) {cal_cluster_e[layer] = e;}
  float get_cal_cluster_e(int layer) const {return cal_cluster_e[layer];}

  float get_x() const{return x;}
  void set_x(float val){x = val;}
  float get_y() const{return y;}
  void set_y(float val){y = val;}
  float get_z() const{return z;}
  void set_z(float val){z = val;}

 private: 

  int     _track_id;
  float   _phi,_d,_kappa,_z0,_dzdl;
  int     clusterID[100];
  float   position[100][3];
  float   momentum;
  float   mom3[3];
  int     charge;
  bool    isprimary;
  bool    ispositive;
  float   quality;
  float   chisq;
  float   chisqv;
  int     ndf;
  float   DCA;
  float   DCA2D;
  float   DCA2Dsigma;
  float   scatter[100];
  float   x,y,z;
  
  TMatrix covariance;
  
  // calorimeter matches
  float   cal_dphi[4];
  float   cal_deta[4];
  float   cal_energy_3x3[4];
  int     cal_cluster_id[4];
  float   cal_cluster_e[4];

  // cluster ids
  //std::multimap<unsigned int,unsigned int> _cluster_ids; //< layer index => cluster id
  
  ClassDef(SvtxTrack,1)
};

#endif

