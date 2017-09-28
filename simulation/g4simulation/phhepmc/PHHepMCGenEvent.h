#ifndef __PHHEPMCGENEVENT__
#define __PHHEPMCGENEVENT__

#include <phool/phool.h>
#include <phool/PHObject.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>


namespace HepMC
{
  class GenEvent;
};


class PHHepMCGenEvent : public PHObject
{
 public:
  
  PHHepMCGenEvent(const int theMomentum = HepMC::Units::GEV,
		  const int theDistance = HepMC::Units::CM);
  PHHepMCGenEvent(const PHHepMCGenEvent& event);
  PHHepMCGenEvent& operator=(const PHHepMCGenEvent& event);
  virtual ~PHHepMCGenEvent();

  virtual void identify(std::ostream& os=std::cout) const;
  virtual void Reset();
  virtual int isValid() const { PHOOL_VIRTUAL_WARNING; return 0; }
  PHHepMCGenEvent* Clone() const {return new PHHepMCGenEvent(*this);}
  
  virtual HepMC::GenEvent* getEvent();
  virtual const HepMC::GenEvent* getEvent() const;

  unsigned int get_id() const {return _id;}
  void set_id(unsigned int id) {_id = id;}
  
  bool addEvent(HepMC::GenEvent *evt);
  bool addEvent(HepMC::GenEvent &evt);
  bool swapEvent(HepMC::GenEvent *evt);
  void clearEvent();

  virtual void moveVertex(double x, double y, double z, double t = 0);

  // the number of entries in the array of particles
  virtual int size(void) const ;
  virtual int vertexSize(void) const ;

  virtual void print(std::ostream& os=std::cout) const;

  int get_momentumunit() const {return _theMomentumUnit;}
  int get_lengthunit() const {return _theDistanceUnit;}

  void PrintEvent();

  bool is_shift_applied() const {return _isVtxShiftApplied;}

protected:

  unsigned int _id;
  bool _isVtxShiftApplied;
  int _theMomentumUnit;
  int _theDistanceUnit;
  HepMC::GenEvent *_theEvt;
  
private:
  

  ClassDef(PHHepMCGenEvent,1)
    
};

#endif	// __PHHEPMCEVENT__
