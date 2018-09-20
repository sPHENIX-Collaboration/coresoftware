#include "PHGenEventv1.h"

#include <HepMC/GenEvent.h>

#include <sstream>

//ClassImp(PHGenEventv1)

using namespace std;

PHGenEventv1::PHGenEventv1()
  : _id(0),
    _event_record(),
    _stale(true),
    _event(NULL)
{}

PHGenEventv1::PHGenEventv1(const unsigned int id, HepMC::GenEvent& event)
   : _id(id),
     _event_record(),
     _stale(true),
     _event(NULL) {  
   set_event(event);
 }

PHGenEventv1::PHGenEventv1(const PHGenEventv1& phevent)
   : _id(phevent.get_id()),
     _event_record(phevent.get_event_record()),
     _stale(true),
     _event(NULL)
{}

PHGenEventv1::PHGenEventv1(const PHGenEventv1* phevent)
   : _id(phevent->get_id()),
     _event_record(phevent->get_event_record()),
     _stale(true),
     _event(NULL)
{}

PHGenEventv1::~PHGenEventv1() {
  if (_event) delete _event;  
}

const HepMC::GenEvent* PHGenEventv1::get_event() const {
  if (stale()) refresh();
  return (const HepMC::GenEvent*)_event;
}

HepMC::GenEvent* PHGenEventv1::get_event() {
  if (stale()) refresh();
  return _event;
}

void PHGenEventv1::set_event(HepMC::GenEvent& event) {

  _event_record.Clear();
  if (_event) {
     delete _event;
     _event = NULL;
  }

  std::stringstream streamer;
//  event.write(streamer);
  _event_record = streamer.str();
  
  refresh();
}

void PHGenEventv1::set_event(HepMC::GenEvent* event) {

  _event_record.Clear();
  if (_event) {
    delete _event;
    _event = NULL;
  }

  std::stringstream streamer;
//  event->write(streamer);
  _event_record = streamer.str();
    
  refresh();
}

size_t PHGenEventv1::particles_size() const {
  if (stale()) refresh();
  return _event->particles_size();
}

size_t PHGenEventv1::vertices_size() const {
  if (stale()) refresh();
  return _event->vertices_size();
}

void PHGenEventv1::Reset() {
  _id = 0;
  _event_record.Clear();
  _stale = true;
  if (_event) {
    delete _event;
    _event = NULL;
  }
}

void PHGenEventv1::print(std::ostream& out) const {
  if (stale()) refresh();
  identify(out);
  out << " id = " << _id << endl;
//  _event->print(out);
}

void PHGenEventv1::refresh() const {

  if (_event) {
    delete _event;
    _event = NULL;
  }

  _event = new HepMC::GenEvent();
  
  std::stringstream streamer;
  streamer << _event_record;
//  _event->read(streamer);
  _stale = false;
}
