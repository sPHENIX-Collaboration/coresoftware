#include <PHHepMCGenEvent.h>
#include "HepMC/GenEvent.h"

#include <TBuffer.h>
#include <TClass.h>

#include <RVersion.h> // root version 

#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <boost/foreach.hpp>

ClassImp(PHHepMCGenEvent)

using namespace std;

PHHepMCGenEvent::PHHepMCGenEvent(const int theMomentum,const int theDistance):
  _theEvt(NULL),
  _isVtxShiftApplied(false),
  _theMomentumUnit(theMomentum),
  _theDistanceUnit(theDistance)
{}


PHHepMCGenEvent::~PHHepMCGenEvent()
{
  if(_theEvt) delete _theEvt;
}


HepMC::GenEvent* PHHepMCGenEvent::getEvent()
{
  return _theEvt;
}

bool PHHepMCGenEvent::addEvent(HepMC::GenEvent *evt)
{
  _theEvt = evt;
  if(!_theEvt) return false;
  return true;
}

bool PHHepMCGenEvent::swapEvent(HepMC::GenEvent *evt)
{
  //if(_theEvt) _theEvt = NULL; 
  _theEvt = evt;

  if(!_theEvt) return false;
  return true;
}


bool PHHepMCGenEvent::addEvent(HepMC::GenEvent &evt)
{
  _theEvt->clear();
  HepMC::GenEvent tmp(evt);
  _theEvt->swap(tmp);
  if(!_theEvt) return false;
  return true;
}

void PHHepMCGenEvent::clearEvent()
{
  if(_theEvt) _theEvt->clear();
}


void PHHepMCGenEvent::moveVertex(double x,double y,double z,double t)
{

  static const float CM2MM = 10.;

  if(!_isVtxShiftApplied)
    {
      for ( HepMC::GenEvent::vertex_iterator vt = _theEvt->vertices_begin();
	    vt != _theEvt->vertices_end(); ++vt )
	{
	  double xShift = (*vt)->position().x() + x*CM2MM;
	  double yShift = (*vt)->position().y() + y*CM2MM;
	  double zShift = (*vt)->position().z() + z*CM2MM;
	  double tShift = (*vt)->position().t() + t;
	  //std::cout << " vertex (x,y,z)= " << x <<" " << y << " " << z << std::endl;
	  (*vt)->set_position( HepMC::FourVector(xShift,yShift,zShift,tShift) ) ;      
	}
      
      _isVtxShiftApplied = true;
    }
  else{ cout << "PHHepMCGenEvent::moveVertex - vertex has already been shifted for this event!" << endl;}

}

int PHHepMCGenEvent::size(void) const
{ 
  return _theEvt->particles_size();
}

int PHHepMCGenEvent::vertexSize(void) const
{ 
  return _theEvt->vertices_size();
}

void PHHepMCGenEvent::Reset()
{
  _isVtxShiftApplied = false;
}

//_____________________________________________________________________________
void PHHepMCGenEvent::identify(std::ostream& os ) const
{
  os << "identify yourself: PHHepMCGenEvent Object" << endl;
  os << "No of Particles: " << size() << endl;
  os << "No of Vertices:  " << vertexSize() << endl;
  return;
}

void PHHepMCGenEvent::print(std::ostream& out) const
{
  identify(out);
}

void PHHepMCGenEvent::PrintEvent() 
{
  _theEvt->print();
}


