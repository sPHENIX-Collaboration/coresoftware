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
  
  PHHepMCGenEvent(const int theMomentum = HepMC::Units::GEV,const int theDistance = HepMC::Units::CM);
  virtual ~PHHepMCGenEvent();
  
  virtual HepMC::GenEvent* getEvent();

  bool addEvent(HepMC::GenEvent *evt);
  bool addEvent(HepMC::GenEvent &evt);
  bool swapEvent(HepMC::GenEvent *evt);
  void clearEvent();

  virtual void moveVertex(double x, double y, double z, double t = 0);

  // the number of entries in the array of particles
  virtual int size(void) const ;
  virtual int vertexSize(void) const ;

  virtual int isValid() const
  { PHOOL_VIRTUAL_WARNING; return 0; }

  virtual void Reset();

  virtual void identify(std::ostream& os=std::cout) const;

  virtual void print(std::ostream& os=std::cout) const;

  int get_momentumunit() const {return _theMomentumUnit;}
  int get_lengthunit() const {return _theDistanceUnit;}

  void PrintEvent();

protected:
  
  HepMC::GenEvent *_theEvt;
  bool _isVtxShiftApplied;
  int _theMomentumUnit;
  int _theDistanceUnit;

private:
  

  ClassDef(PHHepMCGenEvent,1)
    
};

#endif	// __PHHEPMCEVENT__
