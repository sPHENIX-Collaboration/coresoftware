#include "PHHepMCGenEvent.h"

#include <HepMC/GenEvent.h>
#include <HepMC/SimpleVector.h>  // for FourVector

#include <sstream>
#include <utility>               // for swap

using namespace std;

PHHepMCGenEvent::PHHepMCGenEvent()
  : _embedding_id(0)
  , _isSimulated(false)
  , _collisionVertex(0, 0, 0, 0)
  , _theEvt(nullptr)
{
}

PHHepMCGenEvent::PHHepMCGenEvent(const PHHepMCGenEvent& event)
  : _embedding_id(event.get_embedding_id())
  , _isSimulated(event.is_simulated())
  , _collisionVertex(event.get_collision_vertex())
  , _theEvt(nullptr)
{
  _theEvt = new HepMC::GenEvent(*event.getEvent());
  return;
}

PHHepMCGenEvent& PHHepMCGenEvent::operator=(const PHHepMCGenEvent& event)
{
  if (&event == this) return *this;

  Reset();

  _embedding_id = event.get_embedding_id();
  _isSimulated = event.is_simulated();

  const HepMC::GenEvent* hepmc = event.getEvent();
  _theEvt = new HepMC::GenEvent(*(hepmc));

  return *this;
}

PHHepMCGenEvent::~PHHepMCGenEvent()
{
  Reset();
}

void PHHepMCGenEvent::Reset()
{
  _embedding_id = 0;
  _isSimulated = false;
  _collisionVertex.set(0, 0, 0, 0);
  delete _theEvt;
  _theEvt = nullptr;
}

HepMC::GenEvent* PHHepMCGenEvent::getEvent()
{
  return _theEvt;
}

const HepMC::GenEvent* PHHepMCGenEvent::getEvent() const
{
  return _theEvt;
}

bool PHHepMCGenEvent::addEvent(HepMC::GenEvent* evt)
{
  if (_theEvt) delete _theEvt;

  _theEvt = evt;
  if (!_theEvt) return false;
  return true;
}

bool PHHepMCGenEvent::swapEvent(HepMC::GenEvent*& evt)
{
  swap(_theEvt, evt);

  if (!_theEvt) return false;
  return true;
}

bool PHHepMCGenEvent::addEvent(HepMC::GenEvent& evt)
{
  return addEvent(new HepMC::GenEvent(evt));
}

void PHHepMCGenEvent::clearEvent()
{
  if (_theEvt) _theEvt->clear();
}

void PHHepMCGenEvent::moveVertex(double x, double y, double z, double t)
{
  _collisionVertex.setX(_collisionVertex.x() + x);
  _collisionVertex.setY(_collisionVertex.y() + y);
  _collisionVertex.setZ(_collisionVertex.z() + z);
  _collisionVertex.setT(_collisionVertex.t() + t);
}

int PHHepMCGenEvent::size(void) const
{
  if (_theEvt)
    return _theEvt->particles_size();
  else
    return 0;
}

int PHHepMCGenEvent::vertexSize(void) const
{
  if (_theEvt)
    return _theEvt->vertices_size();
  else
    return 0;
}

//_____________________________________________________________________________
void PHHepMCGenEvent::identify(std::ostream& os) const
{
  os << "identify yourself: PHHepMCGenEvent Object, ";
  os << " embedding_id = " << _embedding_id;
  os << " isSimulated = " << _isSimulated;
  os << " collisionVertex = (" << _collisionVertex.x() << "," << _collisionVertex.y() << "," << _collisionVertex.z() << ") cm, " << _collisionVertex.t() << " ns";
  os << ", No of Particles: " << size();
  os << ", No of Vertices:  " << vertexSize();
  os << endl;
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
