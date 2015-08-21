#ifndef __SVTXTRACK_H__
#define __SVTXTRACK_H__

#include <phool/PHObject.h>
#include <TMatrix.h>
#include <iostream>

#include <map>
#include <stdint.h>
#include <climits>
#include <string>

class SvtxTrack : public PHObject
{
  
 public:

  enum CAL_LAYER {PRES,CEMC,HCALIN,HCALOUT};

  SvtxTrack();
  SvtxTrack(SvtxTrack *track);
  SvtxTrack(const SvtxTrack& track);
  virtual ~SvtxTrack() {};
  
  // The "standard PHObject response" functions...
  void identify(std::ostream &os=std::cout) const;
  void Reset();
  int  isValid() const;

  void set_id(int id) {trackID = id;}
  void setTrackID(int index){trackID = index;}
  int getTrackID() const {return trackID;}  
  
  void setClusterID(int layer, int index);
  int getClusterID(int layer) const;
  bool hasCluster(int layer) const;
  
  void setScatter(int layer, float sct);
  float getScatter(int layer) const;
  
  void setHitPosition(int layer, float x, float y, float z);
  float getHitPosition(int layer, int coor) const;
  
  void setMomentum(float p);
  float getMomentum() const;
  
  void set3Momentum(float px, float py, float pz);
  float get3Momentum(int coor) const;
  
  void setCharge(int c);
  int getCharge() const;
  
  void setPrimary(bool prim);
  bool getPrimary() const;
  
  void setPositive(bool prim);
  bool getPositive() const;
  
  //void setNhits(int layer, short n);
  short getNhits() const;
  
  void setQuality(float q);
  float getQuality() const;

  void setChisq(float q);
  float getChisq() const;

  void setChisqv(float q);
  float getChisqv() const;

  void setNDF(int q);
  int getNDF() const;

  void setDCA(float d);
  float getDCA() const;

  void setDCA2D(float d);
  float getDCA2D() const;
  
  void setDCA2Dsigma(float d);
  float getDCA2Dsigma() const;

  float getInnerMostHitPosition(int coor) const;
  
  float phi,d,kappa,z0,dzdl;
  
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

  float get_x() const{return get_property_float(prop_x);}
  void set_x(float val){set_property(prop_x, val);}
  float get_y() const{return get_property_float(prop_y);}
  void set_y(float val){set_property(prop_y, val);}
  float get_z() const{return get_property_float(prop_z);}
  void set_z(float val){set_property(prop_z, val);}

  

 public:

  //! In order to add a new property, please add a short desc. to SvtxTrack::get_property_info()
  enum PROPERTY
  {//
    //! Truth ID
    prop_FastSim_TruthID = 0,

    //!Track properties
    prop_momentum = 10,
    prop_charge,
    prop_quality,
    prop_chisq,
    prop_chisqv,
    prop_ndf,

    //! vertex properties
    prop_DCA = 30,
    prop_DCA2D,
    prop_DCA2Dsigma,
    prop_x,prop_y,prop_z,

    //! max limit in order to fit into 8 bit unsigned number
    prop_MAX_NUMBER = UCHAR_MAX
  };

  enum PROPERTY_TYPE
  {//
    type_int = 1,
    type_uint = 2,
    type_float = 3,
    type_unknown = -1
  };

  virtual bool  has_property(const PROPERTY prop_id) const ;
  virtual float get_property_float(const PROPERTY prop_id) const;
  virtual int   get_property_int(const PROPERTY prop_id) const;
  virtual unsigned int   get_property_uint(const PROPERTY prop_id) const;
  virtual void  set_property(const PROPERTY prop_id, const float value) ;
  virtual void  set_property(const PROPERTY prop_id, const int value) ;
  virtual void  set_property(const PROPERTY prop_id, const unsigned int value) ;
  static std::pair<const std::string,PROPERTY_TYPE> get_property_info(PROPERTY prop_id);
  static bool check_property(const PROPERTY prop_id, const PROPERTY_TYPE prop_type);
  static std::string get_property_type(const PROPERTY_TYPE prop_type);

 protected:

  virtual unsigned int get_property_nocheck(const PROPERTY prop_id) const ;
  void set_property_nocheck(const PROPERTY prop_id,const unsigned int ui) {prop_map[prop_id]=ui;}

  //! storage types for additional property
  typedef uint8_t prop_id_t;
  typedef uint32_t prop_storage_t;
  typedef std::map<prop_id_t, prop_storage_t> prop_map_t;

  //! convert between 32bit inputs and storage type prop_storage_t
  union u_property{
    float fdata;
    int32_t idata;
    uint32_t uidata;

    u_property(int32_t in): idata(in) {}
    u_property(uint32_t in): uidata(in) {}
    u_property(float in): fdata(in) {}
    u_property(): uidata(0) {}

    prop_storage_t get_storage() const {return uidata;}
  };

  //! container for additional property
  prop_map_t prop_map;

  ///@}



 protected:

  int     clusterID[100];
  int     trackID;
  float   position[100][3];
//  float   momentum;
  float   mom3[3];
//  int     charge;
  bool    isprimary;
  bool    ispositive;
//  float   quality;
//  float   chisq;
//  float   chisqv;
//  int     ndf;
//  float   DCA;
//  float   DCA2D;
//  float   DCA2Dsigma;
  float   scatter[100];
//  float   x,y,z;

  TMatrix covariance;

  // calorimeter matches
  float   cal_dphi[4];
  float   cal_deta[4];
  float   cal_energy_3x3[4];
  int     cal_cluster_id[4];
  float   cal_cluster_e[4];


  /** @name Property Tags
   *  Support a variable size propery tags
   */
  ///@{




  ClassDef(SvtxTrack,2)
};


inline void SvtxTrack::setClusterID(int layer, int index){clusterID[layer]=index;}
inline int SvtxTrack::getClusterID(int layer) const {return clusterID[layer];}
inline bool SvtxTrack::hasCluster(int layer) const {return (clusterID[layer]>-9999);}

inline void SvtxTrack::setScatter(int layer, float sct){scatter[layer]=sct;}
inline float SvtxTrack::getScatter(int layer) const {return scatter[layer];}

inline void SvtxTrack::setDCA(float d){set_property(prop_DCA, d);}
inline float SvtxTrack::getDCA() const {return get_property_float(prop_DCA);}
inline void SvtxTrack::setDCA2D(float d){set_property(prop_DCA2D, d);}
inline float SvtxTrack::getDCA2D() const{return get_property_float(prop_DCA2D);}
inline void SvtxTrack::setDCA2Dsigma(float s){set_property(prop_DCA2Dsigma, s);}
inline float SvtxTrack::getDCA2Dsigma() const {return get_property_float(prop_DCA2Dsigma);}

inline void SvtxTrack::setMomentum(float p){set_property(prop_momentum, p);}
inline float SvtxTrack::getMomentum() const {return get_property_float(prop_momentum);}

inline void SvtxTrack::setQuality(float q){set_property(prop_quality, q);}
inline float SvtxTrack::getQuality() const {return get_property_float(prop_quality);}

inline void SvtxTrack::setChisq(float q){set_property(prop_chisq, q);}
inline float SvtxTrack::getChisq() const{return get_property_float(prop_chisq);}

inline void SvtxTrack::setChisqv(float q){set_property(prop_chisqv, q);}
inline float SvtxTrack::getChisqv() const {return get_property_float(prop_chisqv);}

inline void SvtxTrack::setNDF(int q){set_property(prop_ndf, q);}
inline int SvtxTrack::getNDF() const{return get_property_int(prop_ndf);}


inline void SvtxTrack::set3Momentum(float px, float py, float pz)
{
  mom3[0] = px;
  mom3[1] = py;
  mom3[2] = pz;
}
inline float SvtxTrack::get3Momentum(int coor) const {return mom3[coor];}

inline void SvtxTrack::setCharge(int c){set_property(prop_charge, c);}
inline int SvtxTrack::getCharge() const {return get_property_int(prop_charge);}

inline void SvtxTrack::setPrimary(bool prim){isprimary=prim;}
inline bool SvtxTrack::getPrimary() const {return isprimary;}

inline void SvtxTrack::setPositive(bool p){ispositive=p;}
inline bool SvtxTrack::getPositive() const {return ispositive;}

inline float SvtxTrack::getHitPosition(int layer, int coor) const {return position[layer][coor];}
inline void SvtxTrack::setHitPosition(int layer, float x, float y, float z)
{
  position[layer][0]=x;
  position[layer][1]=y;
  position[layer][2]=z;
}


#endif

