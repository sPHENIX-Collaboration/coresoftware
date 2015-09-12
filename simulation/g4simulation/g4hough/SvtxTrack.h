#ifndef __SVTXTRACK_H__
#define __SVTXTRACK_H__

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <cmath>

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

  void setMomentum(float p) {}
  float getMomentum() const {
    float px = _states.find(0.0)->second.get_px();
    float py = _states.find(0.0)->second.get_py();
    float pz = _states.find(0.0)->second.get_pz();
    return sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
  }
  
  void set3Momentum(float px, float py, float pz) {
    _states[0.0].set_px(px);
    _states[0.0].set_py(py);
    _states[0.0].set_pz(pz);
  };
  float get3Momentum(int coor) const {
    return _states.find(0.0)->second.get_mom(coor);
  }
  
  void setCharge(int c) {
    if (c > 0) setPositive(true);
    else setPositive(false);
  }
  int getCharge() const {
    if (getPositive()) return +1;
    else return -1;
  }
  
  void setPositive(bool prim) {_is_positive_charge = prim;}
  bool getPositive() const {return _is_positive_charge;}

  void setPrimary(bool prim) {
    std::cout << "SvtxTrack:: ERROR - deprecated interface call" << std::endl;
  }
  bool getPrimary() const {
    return false;
  }
  
  //void setNhits(int layer, short n);
  short getNhits() const;
  
  void setQuality(float q) {}
  float getQuality() const {
    if (_ndf!=0) return _chisq/_ndf;
    return NAN;
  }

  void setChisq(float q) {_chisq = q;}
  float getChisq() const {return _chisq;}

  void setChisqv(float q) {
    std::cout << "SvtxTrack:: ERROR - deprecated interface call" << std::endl;
  }
  float getChisqv() const {
    std::cout << "SvtxTrack:: ERROR - deprecated interface call" << std::endl;
    return NAN;
  }

  void setNDF(int q) {_ndf = q;}
  int  getNDF() const {return _ndf;}

  void  setDCA(float d) {_DCA = d;}
  float getDCA() const {return _DCA;}

  void  setDCA2D(float d) {_DCA2D = d;}
  float getDCA2D() const {return _DCA2D;}
  
  void  setDCA2Dsigma(float s) {_DCA2Dsigma = s;}
  float getDCA2Dsigma() const {return _DCA2Dsigma;}

  float getInnerMostHitPosition(int coor) const;

  void  set_phi(float phi) {_states[0.0].set_phi(phi);}
  float get_phi() const {return _states.find(0.0)->second.get_phi();}
    
  void  set_d(float d) {_states[0.0].set_d(d);}
  float get_d() const {return _states.find(0.0)->second.get_d();}
  
  void  set_kappa(float kappa) {_states[0.0].set_kappa(kappa);}
  float get_kappa() const {return _states.find(0.0)->second.get_kappa();}
    
  void  set_z0(float z0) {_states[0.0].set_z0(z0);}
  float get_z0() const {return _states.find(0.0)->second.get_z0();}
    
  void  set_dzdl(float dzdl) {_states[0.0].set_dzdl(dzdl);}
  float get_dzdl() const {return _states.find(0.0)->second.get_dzdl();}
  
  float getCovariance(int i,int j) const {return get_error(i,j);}
  void  setCovariance(int i,int j, float val) {set_error(i,j,val);}

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

  float get_x() const  {return _states.find(0.0)->second.get_x();}
  void  set_x(float x) {_states[0.0].set_x(x);}
  
  float get_y() const  {return _states.find(0.0)->second.get_y();}
  void  set_y(float y) {_states[0.0].set_y(y);}

  float get_z() const  {return _states.find(0.0)->second.get_z();}
  void  set_z(float z) {_states[0.0].set_z(z);}

  // add convience calculations
  //float get_eta() const;
  //float get_theta() const;
  //float get_phi() const;
  //float get_pt() const;
  //float get_p() const;

 private: 

  // keep it private for now to minimize interface changes
  // --- inner State class ---------------------------------------------------//
  class State {                                                               //
  public:                                                                     //
    State();                                                                  //
    virtual ~State() {}                                                       //
                                                                              //
    float get_x() const {return _x;}                                          //
    void  set_x(float x) {_x = x;}                                            //
                                                                              //
    float get_y() const {return _y;}                                          //
    void  set_y(float y) {_y = y;}                                            //
                                                                              //
    float get_z() const {return _z;}                                          //
    void  set_z(float z) {_z = z;}                                            //
                                                                              //
    float get_px() const {return _mom[0];}                                    //
    void  set_px(float px) {_mom[0] = px;}                                    //
                                                                              //
    float get_py() const {return _mom[1];}                                    //
    void  set_py(float py) {_mom[1] = py;}                                    //
                                                                              //
    float get_pz() const {return _mom[2];}                                    //
    void  set_pz(float pz) {_mom[2] = pz;}                                    //
                                                                              //
    float get_mom(unsigned int i) const {return _mom[i];}                     //
                                                                              //
    float get_error(int i, int j) const;                                      //
    void  set_error(int i, int j, float value);                               //
                                                                              //
    void  set_phi(float phi) {_phi = phi;}                                    //
    float get_phi() const {return _phi;}                                      //
                                                                              //
    void  set_d(float d) {_d = d;}                                            //
    float get_d() const {return _d;}                                          //
                                                                              //
    void  set_kappa(float kappa) {_kappa = kappa;}                            //
    float get_kappa() const {return _kappa;}                                  //
                                                                              //
    void  set_z0(float z0) {_z0 = z0;}                                        //
    float get_z0() const {return _z0;}                                        //
                                                                              //
    void  set_dzdl(float dzdl) {_dzdl = dzdl;}                                //
    float get_dzdl() const {return _dzdl;}                                    //
                                                                              //
  private:                                                                    //
    float _x;                                                                 //
    float _y;                                                                 //
    float _z;                                                                 //
    float _mom[3];                                                            //
    std::vector<std::vector<float> > _covar;                                  //
    float _phi;                                                               //
    float _d;                                                                 //
    float _kappa;                                                             //
    float _z0;                                                                //
    float _dzdl;                                                              //
  };                                                                          //
  // --- inner State class ---------------------------------------------------//

  // class CaloMatch
  
  // keep these private for now
  // attempting ~zero interface changes during refactor
  float get_error(int i, int j) const {return _states.find(0.0)->second.get_error(i,j);}
  void  set_error(int i, int j, float value) {return _states[0.0].set_error(i,j,value);}
  
  // track information
  unsigned int _track_id;
  bool         _is_positive_charge;
  float        _chisq;
  unsigned int _ndf;

  // extended track information (non-primary tracks only)
  float _DCA;
  float _DCA2D;
  float _DCA2Dsigma;

  // extended track information (primary tracks only)
  // unsigned int _vertex_id;
  
  // track state information
  std::map<float,SvtxTrack::State> _states; //< path length => state object
  
  // cluster contents
  std::map<int,unsigned int> _cluster_ids; //< layer index => cluster id

  // the cluster positions aren't really useful on their own
  // without the cluster uncertainties... maybe we should eliminate
  // this member for further storage gains (and use the ids to fetch the clusters
  // for remaking the fits)
  // we will first need to replace the public projection method to use and outer
  // state vector instead of the hit postion
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

