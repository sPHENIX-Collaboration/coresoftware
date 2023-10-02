// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTPLANEINFO_H
#define EVENTPLANEINFO_H

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>
#include <map>

class Eventplaneinfo : public PHObject
{
 public:
  
  enum EPTYPE
  {
    UNDEFINED = 0,
    EPDS = 10,
    EPDN = 20,
    MBDS = 30,
    MBDN = 40
  };

  typedef std::map<Eventplaneinfo::EPTYPE, unsigned int>::const_iterator ConstEpIter;
  typedef std::map<Eventplaneinfo::EPTYPE, unsigned int>::iterator EpIter;

  ~Eventplaneinfo() override {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Eventplaneinfo base class" << std::endl;
  }
  PHObject* CloneMe() const override { return nullptr; }

     
  virtual unsigned int get_id() const { return 0xFFFFFFFF; }
  virtual void set_id(unsigned int) {}
    
  virtual double get_psi_raw(int /*order*/) const {return NAN;}
  virtual void set_psi_raw(unsigned int /*order*/, double /*f*/) {return;}

  virtual bool empty_ep_ids() const { return true; }
  virtual size_t size_ep_ids() const { return 0; }
  virtual size_t count_ep_ids(EPTYPE /*type*/) const { return 0; }

  virtual void clear_ep_ids() {}
  virtual void insert_ep_ids(EPTYPE /*type*/, unsigned int /*ep_id*/) {}
  virtual size_t erase_ep_ids(EPTYPE /*type*/) { return 0; }
  virtual void erase_ep_ids(EpIter /*iter*/) {}
  virtual void erase_ep_ids(EpIter /*first*/, EpIter /*last*/) {}

  virtual ConstEpIter begin_ep_ids() const;
  virtual ConstEpIter find_ep_ids(EPTYPE type) const;
  virtual ConstEpIter end_ep_ids() const;

  virtual EpIter begin_ep_ids();
  virtual EpIter find_ep_ids(EPTYPE type);
  virtual EpIter end_ep_ids();

 protected:
  Eventplaneinfo() {}

 private:
  ClassDefOverride(Eventplaneinfo, 1);
};

#endif
