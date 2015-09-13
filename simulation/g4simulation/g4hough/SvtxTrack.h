#ifndef __SVTXTRACK_H__
#define __SVTXTRACK_H__

#include <phool/PHObject.h>

#include <iostream>
#include <set>
#include <map>
#include <cmath>

class SvtxTrack : public PHObject {
  
 public:

  // --- inner State class ---------------------------------------------------//
  class State {                                                               //
  public:                                                                     //
    State(float pathlength = 0.0);                                            //
    virtual ~State() {}                                                       //
                                                                              //
    float get_pathlength() const {return _pathlength;}                        //
                                                                              //
    float get_x() const {return _pos[0];}                                     //
    void  set_x(float x) {_pos[0] = x;}                                       //
                                                                              //
    float get_y() const {return _pos[1];}                                     //
    void  set_y(float y) {_pos[1] = y;}                                       //
                                                                              //
    float get_z() const {return _pos[2];}                                     //
    void  set_z(float z) {_pos[2] = z;}                                       //
                                                                              //
    float get_pos(unsigned int i) const {return _pos[i];}                     //
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
    float get_error(unsigned int i, unsigned int j) const;                    //
    void  set_error(unsigned int i, unsigned int j, float value);             //
                                                                              //
    float get_helix_phi() const {return _helix_phi;}                          //
    void  set_helix_phi(float helix_phi) {_helix_phi = helix_phi;}            //
                                                                              //
    float get_helix_d() const {return _helix_d;}                              //
    void  set_helix_d(float d) {_helix_d = d;}                                //
                                                                              //
    float get_helix_kappa() const {return _helix_kappa;}                      //
    void  set_helix_kappa(float kappa) {_helix_kappa = kappa;}                //
                                                                              //
    float get_helix_z0() const {return _helix_z0;}                            //
    void  set_helix_z0(float z0) {_helix_z0 = z0;}                            //
                                                                              //
    float get_helix_dzdl() const {return _helix_dzdl;}                        //
    void  set_helix_dzdl(float dzdl) {_helix_dzdl = dzdl;}                    //
                                                                              //
  private:                                                                    //
                                                                              //
    unsigned int covar_index(unsigned int i, unsigned int j) const;           //
                                                                              //
    float _pathlength;                                                        //
    float _pos[3];                                                            //
    float _mom[3];                                                            //
    float _covar[21]; // 6x6 triangular packed storage                        //
    float _helix_phi;                                                         //
    float _helix_d;                                                           //
    float _helix_kappa;                                                       //
    float _helix_z0;                                                          //
    float _helix_dzdl;                                                        //
  };                                                                          //
  // --- inner State class ---------------------------------------------------//
 
  typedef std::map<float,SvtxTrack::State>::const_iterator ConstStateIter;
  typedef std::map<float,SvtxTrack::State>::iterator       StateIter; 
  
  typedef std::set<unsigned int>::const_iterator         ConstClusterIter;
  typedef std::set<unsigned int>::iterator               ClusterIter;
  
  enum CAL_LAYER {PRES=0,CEMC=1,HCALIN=2,HCALOUT=3};

  SvtxTrack();
  virtual ~SvtxTrack() {}
  
  // The "standard PHObject response" functions...
  void identify(std::ostream &os=std::cout) const;
  void Reset();
  int  isValid() const;

  unsigned int get_id() const          {return _track_id;}
  void         set_id(unsigned int id) {_track_id = id;}

  bool get_positive_charge() const {return _is_positive_charge;}
  void set_positive_charge(bool ispos) {_is_positive_charge = ispos;}

  int  get_charge() const {return (get_positive_charge()) ? 1 : -1;}
  void set_charge(int charge) {(charge > 0) ? set_positive_charge(true) : set_positive_charge(false);}

  float get_chisq() const {return _chisq;}  
  void  set_chisq(float chisq) {_chisq = chisq;}

  unsigned int get_ndf() const {return _ndf;}
  void         set_ndf(int ndf) {_ndf = ndf;}

  float get_quality() const {return (_ndf != 0) ? _chisq/_ndf : NAN;}

  float get_dca() const {return _dca;}
  void  set_dca(float dca) {_dca = dca;}

  float get_dca2d() const {return _dca2d;}  
  void  set_dca2d(float dca2d) {_dca2d = dca2d;}

  float get_dca2d_error() const {return _dca2d_error;}  
  void  set_dca2d_error(float error) {_dca2d_error = error;}

  float get_x() const  {return _states.find(0.0)->second.get_x();}
  void  set_x(float x) {_states[0.0].set_x(x);}
  
  float get_y() const  {return _states.find(0.0)->second.get_y();}
  void  set_y(float y) {_states[0.0].set_y(y);}

  float get_z() const  {return _states.find(0.0)->second.get_z();}
  void  set_z(float z) {_states[0.0].set_z(z);}

  float get_pos(unsigned int i) const {return _states.find(0.0)->second.get_pos(i);}

  float get_px() const   {return _states.find(0.0)->second.get_px();}
  void  set_px(float px) {_states[0.0].set_px(px);}
  
  float get_py() const   {return _states.find(0.0)->second.get_py();}
  void  set_py(float py) {_states[0.0].set_py(py);}

  float get_pz() const   {return _states.find(0.0)->second.get_pz();}
  void  set_pz(float pz) {_states[0.0].set_pz(pz);}

  float get_mom(unsigned int i) const {return _states.find(0.0)->second.get_mom(i);}

  float get_p() const   {return sqrt(pow(get_px(),2) + pow(get_py(),2) + pow(get_pz(),2));}
  float get_pt() const  {return sqrt(pow(get_px(),2) + pow(get_py(),2));}
  float get_eta() const {return asinh(get_pz()/get_pt());}
  float get_phi() const {return atan2(get_py(),get_px());}

  float get_error(int i, int j) const {return _states.find(0.0)->second.get_error(i,j);}
  void  set_error(int i, int j, float value) {return _states[0.0].set_error(i,j,value);}

  float get_helix_phi() const {return _states.find(0.0)->second.get_helix_phi();}  
  void  set_helix_phi(float phi) {_states[0.0].set_helix_phi(phi);}

  float get_helix_d() const {return _states.find(0.0)->second.get_helix_d();}
  void  set_helix_d(float d) {_states[0.0].set_helix_d(d);}

  float get_helix_kappa() const {return _states.find(0.0)->second.get_helix_kappa();}  
  void  set_helix_kappa(float kappa) {_states[0.0].set_helix_kappa(kappa);}

  float get_helix_z0() const {return _states.find(0.0)->second.get_helix_z0();}    
  void  set_helix_z0(float z0) {_states[0.0].set_helix_z0(z0);}

  float get_helix_dzdl() const {return _states.find(0.0)->second.get_helix_dzdl();}    
  void  set_helix_dzdl(float dzdl) {_states[0.0].set_helix_dzdl(dzdl);}

  //
  // state methods
  //
  bool   empty_states()                 const {return _states.empty();}
  size_t size_states()                  const {return _states.size();}
  size_t count_states(float pathlength) const {return _states.count(pathlength);}
  void   clear_states() {
    _states.clear();
    insert_state(State(0.0));
  }
  
  const State* get_state(float pathlength) const;
        State* get_state(float pathlength); 
        State* insert_state(const State &state);
        size_t erase_state(float pathlength) {
	  _states.erase(pathlength);
	  if (pathlength == 0) insert_state(State(0.0));
	  return _states.size();
	}

  ConstStateIter begin_states()                const {return _states.begin();}
  ConstStateIter  find_state(float pathlength) const {return _states.find(pathlength);}
  ConstStateIter   end_states()                const {return _states.end();}

  StateIter begin_states()                {return _states.begin();}
  StateIter  find_state(float pathlength) {return _states.find(pathlength);}
  StateIter   end_states()                {return _states.end();}
    
  //
  // associated cluster ids methods
  //
  void                clear_clusters()                           {_cluster_ids.clear();}
  bool                empty_clusters() const                     {return _cluster_ids.empty();}
  size_t              size_clusters() const                      {return _cluster_ids.size();}
  void                insert_cluster(unsigned int clusterid)     {_cluster_ids.insert(clusterid);}
  size_t              erase_cluster(unsigned int clusterid)      {return _cluster_ids.erase(clusterid);}
  ConstClusterIter    begin_clusters() const                     {return _cluster_ids.begin();}
  ConstClusterIter    find_cluster(unsigned int clusterid) const {return _cluster_ids.find(clusterid);}
  ConstClusterIter    end_clusters() const                       {return _cluster_ids.end();}
  ClusterIter         begin_clusters()                           {return _cluster_ids.begin();}
  ClusterIter         find_cluster(unsigned int clusterid)       {return _cluster_ids.find(clusterid);}
  ClusterIter         end_clusters()                             {return _cluster_ids.end();}

  //
  // calo projection methods
  //
  void  set_cal_energy_3x3(CAL_LAYER layer, float energy_3x3);
  float get_cal_energy_3x3(CAL_LAYER layer) const;

  void         set_cal_cluster_id(CAL_LAYER layer, unsigned int id);
  unsigned int get_cal_cluster_id(CAL_LAYER layer) const;
  
  void  set_cal_dphi(CAL_LAYER layer, float dphi);
  float get_cal_dphi(CAL_LAYER layer) const;

  void  set_cal_deta(CAL_LAYER layer, float deta);
  float get_cal_deta(CAL_LAYER layer) const;

  void  set_cal_cluster_e(CAL_LAYER layer, float e);
  float get_cal_cluster_e(CAL_LAYER layer) const;
  
 private: 

  // keep it private for now to minimize interface changes
  // --- inner CaloProjection class ------------------------------------------//
  class CaloProjection {                                                      //
  public:                                                                     //
    CaloProjection();                                                         //
    virtual ~CaloProjection() {}                                              //
                                                                              //
    float get_energy_3x3() const {return _e3x3;}                              //
    void  set_energy_3x3(float e3x3) {_e3x3 = e3x3;}                          //
                                                                              //
    unsigned int get_cluster_id() const {return _clus_id;}                    //
    void         set_cluster_id(unsigned int clus_id) {_clus_id = clus_id;}   //
                                                                              //
    float get_deta() const {return _deta;}                                    //
    void  set_deta(float deta) {_deta = deta;}                                //
                                                                              //
    float get_dphi() const {return _dphi;}                                    //
    void  set_dphi(float dphi) {_dphi = dphi;}                                //
                                                                              //
    float get_cluster_energy() const {return _clus_e;}                        //
    void  set_cluster_energy(float clus_e) {_clus_e = clus_e;}                //
                                                                              //
  private:                                                                    //
    float _e3x3;                                                              //
    unsigned int _clus_id;                                                    //
    float _deta;                                                              //
    float _dphi;                                                              //
    float _clus_e;                                                            //
  };                                                                          //
  // --- inner CaloProjection class ------------------------------------------//
  
  // track information
  unsigned int _track_id;
  bool         _is_positive_charge;
  float        _chisq;
  unsigned int _ndf;

  // extended track information (non-primary tracks only)
  float _dca;
  float _dca2d;
  float _dca2d_error;

  // extended track information (primary tracks only)
  // unsigned int _vertex_id;
  
  // track state information
  std::map<float,SvtxTrack::State> _states; //< path length => state object
  
  // cluster contents
  std::set<unsigned int> _cluster_ids;
  
  // calorimeter matches
  std::map<CAL_LAYER,SvtxTrack::CaloProjection> _calo_matches;
  
  ClassDef(SvtxTrack,1)
};

#endif

