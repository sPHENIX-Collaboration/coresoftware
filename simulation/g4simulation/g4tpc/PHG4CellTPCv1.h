#ifndef PHG4CellTPCv1_h__
#define PHG4CellTPCv1_h__

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellDefs.h>
#ifdef __CINT__
#include <stdint.h>
#else
#include <cstdint>
#endif
#include <iostream>
#include <map>

class PHG4CellTPCv1: public PHG4Cell
{
 public:
  PHG4CellTPCv1();
  PHG4CellTPCv1(const PHG4CellDefs::keytype g4cellid);
  virtual ~PHG4CellTPCv1();

  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset();

  void set_cellid(const PHG4CellDefs::keytype i) {cellid = i;}

  PHG4CellDefs::keytype get_cellid() const {return cellid;}
  bool has_binning(const PHG4CellDefs::CellBinning binning) const;
  short int get_detid() const;

  void add_edep(const PHG4HitDefs::keytype g4hitid, const int timebin, const float edep);

  /* EdepConstRange get_g4hits() { */
  /*   return std::make_pair(hitedeps.begin(), hitedeps.end()); */
  /* } */
  


  void print() const;

  bool  has_property(const PROPERTY prop_id) const;
  float get_property_float(const PROPERTY prop_id) const;
  int   get_property_int(const PROPERTY prop_id) const;
  unsigned int   get_property_uint(const PROPERTY prop_id) const;
  void  add_property(const PROPERTY prop_id, const float value);
  void  add_property(const PROPERTY prop_id, const int value);
  void  add_property(const PROPERTY prop_id, const unsigned int value);
  void  set_property(const PROPERTY prop_id, const float value);
  void  set_property(const PROPERTY prop_id, const int value);
  void  set_property(const PROPERTY prop_id, const unsigned int value);


 protected:
  unsigned int get_property_nocheck(const PROPERTY prop_id) const;
  void set_property_nocheck(const PROPERTY prop_id,const unsigned int ui) {prop_map[prop_id]=ui;}

  PHG4CellDefs::keytype cellid;
  std::map<int,EdepMap> timeseq;
//  tpctod trainOfDigits;

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


  ClassDef(PHG4CellTPCv1,1)
};

#endif
