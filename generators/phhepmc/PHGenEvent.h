#ifndef PHHEPMC_PHGENEVENT_H
#define PHHEPMC_PHGENEVENT_H

#include <phool/PHObject.h>

#include <HepMC/Units.h>
#include <HepMC/GenEvent.h>

class PHGenEvent : public PHObject {
  
 public:
  virtual ~PHGenEvent() {}

  virtual const HepMC::GenEvent* get_event() const {return NULL;}
  virtual HepMC::GenEvent* get_event() {return NULL;}
  virtual void set_event(HepMC::GenEvent& event) {}
  virtual void set_event(HepMC::GenEvent* event) {}

  virtual unsigned int get_id() const {return 0;}
  virtual void set_id(const unsigned int id) {}
  
  // the number of entries in the array of particles
  virtual size_t particles_size() const {return 0;}
  virtual size_t vertices_size() const {return 0;}

  virtual void identify(std::ostream& out = std::cout) const {
    out << "PHGenEvent" << std::endl;
  }
  virtual void Reset() {}
  virtual void print(std::ostream& out = std::cout) const {}

protected:
  PHGenEvent() {}

private:

  ClassDef(PHGenEvent,1)    
};

#endif	// PHHEPMC_PHHEPMCEVENT_H
