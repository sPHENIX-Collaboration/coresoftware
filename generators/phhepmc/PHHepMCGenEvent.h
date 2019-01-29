#ifndef PHHEPMC_PHHEPMCGENEVENT_H
#define PHHEPMC_PHHEPMCGENEVENT_H

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>
#include <HepMC/SimpleVector.h>

namespace HepMC
{
class GenEvent;
};

class PHHepMCGenEvent : public PHObject
{
 public:
  PHHepMCGenEvent();

  PHHepMCGenEvent(const PHHepMCGenEvent& event);
  PHHepMCGenEvent& operator=(const PHHepMCGenEvent& event);
  virtual ~PHHepMCGenEvent();

  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset();
  virtual int isValid() const
  {
    PHOOL_VIRTUAL_WARNING;
    return 0;
  }
  PHHepMCGenEvent* Clone() const { return new PHHepMCGenEvent(*this); }
  virtual HepMC::GenEvent* getEvent();
  virtual const HepMC::GenEvent* getEvent() const;

  //! embedding ID for the event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int get_embedding_id() const { return _embedding_id; }

  //! embedding ID for the event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  void set_embedding_id(int id) { _embedding_id = id; }

  //! whether this event has been processed in Geant4 simulation
  bool is_simulated() const { return _isSimulated; }

  //! whether this event has been processed in Geant4 simulation
  void is_simulated(bool v) { _isSimulated = v; }

  //! collision vertex position in the Hall coordinate system, use PHENIX units of cm, ns
  const HepMC::FourVector& get_collision_vertex() const { return _collisionVertex; }

  //! collision vertex position in the Hall coordinate system, use PHENIX units of cm, ns
  void set_collision_vertex(const HepMC::FourVector& v) { _collisionVertex = v; }

  //! host an HepMC event
  bool addEvent(HepMC::GenEvent* evt);
  bool addEvent(HepMC::GenEvent& evt);
  bool swapEvent(HepMC::GenEvent*& evt);
  void clearEvent();

  //! move the collision vertex position in the Hall coordinate system, use PHENIX units of cm, ns
  virtual void moveVertex(double x, double y, double z, double t = 0);

  // the number of entries in the array of particles
  virtual int size(void) const;
  virtual int vertexSize(void) const;

  virtual void print(std::ostream& os = std::cout) const;

  void PrintEvent();

 protected:
  //! \brief Embedding ID for this generated event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int _embedding_id;

  //! whether this event has been processed in Geant4 simulation
  bool _isSimulated;

  //! collision vertex position in the Hall coordinate system, use PHENIX units of cm, ns
  HepMC::FourVector _collisionVertex;

  //! The HEP MC record from event generator. Note the units are recorded in GenEvent
  HepMC::GenEvent* _theEvt;

  ClassDef(PHHepMCGenEvent, 5)
};

#endif  // PHHEPMC_PHHEPMCEVENT_H
