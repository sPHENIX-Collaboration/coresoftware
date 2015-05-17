#ifndef __JET_H__
#define __JET_H__

#include <phool/PHObject.h>
#include <map>
#include <iostream>

class Jet : public PHObject
{

public:

  /*! \addtogroup Generic Features
   *  @{
   */

  Jet();
  virtual
  ~Jet()
  {
  }

  // PHObject virtual overloads

  virtual void
  identify(std::ostream& os = std::cout) const;
  virtual void
  Reset();
  virtual int
  isValid() const;

  /*! @} */

  // jet info
  virtual unsigned int
  get_id() const;
  virtual void
  set_id(unsigned int id);

  virtual float
  get_px() const;
  virtual void
  set_px(float px);

  virtual float
  get_py() const;
  virtual void
  set_py(float py);

  virtual float
  get_pz() const;
  virtual void
  set_pz(float pz);

  virtual float
  get_e() const;
  virtual void
  set_e(float e);

  /*! @} */

  /*! \addtogroup clustered component
   * clustered component methods (multimap interface based)
   * source type id --> unique id within that storage
   *  @{
   */
  enum ALGO
  {
    NONE, ANTIKT, KT, CONE
  };

  enum SRC
  {
    TRACKS,
    EM_TOWERS,
    EM_CLUSTERS,
    HCALIN_TOWERS,
    HCALIN_CLUSTERS,
    HCALOUT_TOWERS,
    HCALOUT_CLUSTERS
  };

  typedef std::multimap<SRC, unsigned int> typ_comp_ids;
  typedef typ_comp_ids::const_iterator ConstIter;
  typedef typ_comp_ids::iterator Iter;

  virtual bool
  empty_comp() const
  {
    return true;
  }
  virtual size_t
  size_comp() const
  {
    return 0;
  }
  virtual size_t
  count_comp(SRC source) const
  {
    return 0;
  }

  virtual void
  clear_comp()
  {
    return;
  }
  virtual void
  insert_comp(SRC source, unsigned int compid)
  {
    return;
  }
  virtual size_t
  erase_comp(SRC source)
  {
    return 0;
  }
  virtual void
  erase_comp(Iter iter)
  {
    return;
  }
  virtual void
  erase_comp(Iter first, Iter last)
  {
    return;
  }

  virtual ConstIter
  begin_comp() const;

  virtual ConstIter
  lower_bound_comp(SRC source) const;

  virtual ConstIter
  upper_bound_comp(SRC source) const;

  virtual ConstIter
  find(SRC source) const;

  virtual ConstIter
  end_comp() const;

  virtual Iter
  begin_comp();

  virtual Iter
  lower_bound_comp(SRC source);

  virtual Iter
  upper_bound_comp(SRC source);

  virtual Iter
  find(SRC source);

  virtual Iter
  end_comp();

  /*! @} */

  /*! \addtogroup Property Tags
   *  Tag the jet object with various tages
   *
   *  Example to try it out in command lines
   *
   *
<code>
[jinhuang@rcas2067 macros]$ root
root [0] gSystem->Load("libg4jets");
root [1] JetV1 j

root [2] j.identify()
---Jet V1-----------------------
jetid: 4294967295
 (px,py,pz,e) =  (nan, nan, nan, nan) GeV
-----------------------------------------------

root [3] j.set_property(Jet::prop_R , 0.2)
root [5] j.set_property(Jet::prop_BFrac  , 0.5)
root [6] j.identify()
---Jet V1-----------------------
jetid: 4294967295
 (px,py,pz,e) =  (nan, nan, nan, nan) GeV
 Jet Radius = 0.2
 Jet B-quark fraction = 0.5

root [7] j.get_property(Jet::prop_BFrac)
(const float)5.00000000000000000e-01

-----------------------------------------------
</code>
   *
   *
   *  @{
   */

  //! Property ID List
  //! You are welcome to add to this list, but please do not remove or change previous entries
  //! Please add description to JetV1::print_property() for new property tags
  enum PROPERTY
  {

    //! Jet radius
    prop_R = 1,

    //! Jet Mass
    prop_JetMass = 11,

    //! Jet Charge
    prop_JetCharge = 12,

    //! B-jet fraction
    prop_BFrac = 101,

    //! Last property tag
    prop_MaxValue
  };
  /*! @} */


  //! print out all existing properties
  virtual void print_property(ostream& os) const {}

  //! whether a property exists
  virtual
  bool
  has_property(PROPERTY prop_id) const
  {
    return false;
  }

  //! get property value
  virtual
  float
  get_property(PROPERTY prop_id) const;

  //! set property value
  virtual
  void
  set_property(PROPERTY prop_id, float value);

private:
  //! a dummpy comp container to make returns for base-class interface functions
  static typ_comp_ids _dummy_ids;

ClassDef(Jet, 1)
  ;
};

#endif
