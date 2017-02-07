#ifndef PHG4Cell_H__
#define PHG4Cell_H__

#include "PHG4CellDefs.h"
#include <g4main/PHG4Hit.h>
#include <phool/PHObject.h>

#include <cmath>
#include <climits>
#include <map>

class PHG4Cell: public PHObject
{
 public:
  typedef std::map<PHG4HitDefs::keytype, float> EdepMap;
  typedef EdepMap::iterator EdepIterator;
  typedef EdepMap::const_iterator EdepConstIterator;
  typedef std::pair<EdepIterator, EdepIterator> EdepRange;
  typedef std::pair<EdepConstIterator, EdepConstIterator> EdepConstRange;

  typedef std::map<int, float> ShowerEdepMap;
  typedef ShowerEdepMap::iterator ShowerEdepIterator;
  typedef ShowerEdepMap::const_iterator ShowerEdepConstIterator;
  typedef std::pair<ShowerEdepIterator, ShowerEdepIterator> ShowerEdepRange;
  typedef std::pair<ShowerEdepConstIterator, ShowerEdepConstIterator> ShowerEdepConstRange;

  virtual ~PHG4Cell() {}

  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Copy(PHG4Cell const &g4cell);
  friend std::ostream &operator<<(std::ostream & stream, const PHG4Cell * cell);
  virtual void Reset();

  // The indices here represent the entry and exit points of the particle
  virtual unsigned int get_layer() const {return UINT_MAX;}
  virtual double get_edep() const {return NAN;}
  virtual double get_eion() const {return NAN;}
  virtual float get_light_yield() const {return NAN;}
  virtual double get_path_length() const {return NAN;}
  virtual PHG4CellDefs::keytype get_cell_id() const {return UINT_MAX;}

  virtual void set_edep(const float f) {return;}
  virtual void set_eion(const float f) {return;}
  virtual void set_light_yield(const float lightYield){return;}
  virtual void set_layer(const unsigned int i) {return;}
  virtual void set_cell_id(const PHG4CellDefs::keytype i) {return;}

  virtual void print() const {std::cout<<"PHG4Cell base class - print() not implemented"<<std::endl;}


  //! Procedure to add a new PROPERTY tag:
  //! 1.add new tag below with unique value,
  //! 2.add a short name to PHG4Cell::get_property_info
  enum PROPERTY 
  {//
    // first various coordinates 1-20
    //! layer ID
    prop_layer = 1,
    //-- summed energy:  - 20-30  --
    //! deposited energy
    prop_edep = 21,
    //! ionizing energy loss
    prop_eion = 22,

    //! for scintillation detectors, the amount of light produced
    prop_light_yield = 23,

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

  virtual bool  has_property(const PROPERTY prop_id) const {return false;}
  virtual float get_property_float(const PROPERTY prop_id) const {return NAN;}
  virtual int   get_property_int(const PROPERTY prop_id) const {return INT_MIN;}
  virtual unsigned int   get_property_uint(const PROPERTY prop_id) const {return UINT_MAX;}
  virtual void  set_property(const PROPERTY prop_id, const float value) {return;}
  virtual void  set_property(const PROPERTY prop_id, const int value) {return;}
  virtual void  set_property(const PROPERTY prop_id, const unsigned int value) {return;}
  static std::pair<const std::string,PROPERTY_TYPE> get_property_info(PROPERTY prop_id);
  static bool check_property(const PROPERTY prop_id, const PROPERTY_TYPE prop_type);
  static std::string get_property_type(const PROPERTY_TYPE prop_type);

 protected:
  PHG4Cell() {}
  virtual unsigned int get_property_nocheck(const PROPERTY prop_id) const {return UINT_MAX;}
  virtual void set_property_nocheck(const PROPERTY prop_id,const unsigned int) {return;}
  ClassDef(PHG4Cell,1)
};


#endif
