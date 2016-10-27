#include "PHHepMCParticleSelectorDecayProduct.h"

#include "PHHepMCGenEvent.h"

#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>

using namespace std;

PHHepMCParticleSelectorDecayProduct::PHHepMCParticleSelectorDecayProduct(const string &name):
  SubsysReco(name)
{
  _theTrigger = 0;
  _theParticle = 11;
  return;
}

PHHepMCParticleSelectorDecayProduct::~PHHepMCParticleSelectorDecayProduct()
{
  return;
}

int
PHHepMCParticleSelectorDecayProduct::InitRun(PHCompositeNode *topNode)
{
  return 0;
}

int PHHepMCParticleSelectorDecayProduct::process_event(PHCompositeNode *topNode)
{
  PHHepMCGenEvent *inEvent = findNode::getClass<PHHepMCGenEvent>(topNode,"PHHepMCGenEvent");
  if(!inEvent) {cerr << "ERROR: PHHepMCGenEvent node not found!" << endl; return -1;}

  HepMC::GenEvent* event = inEvent->getEvent();
  int npart = event->particles_size();
  int nvert = event->vertices_size();
  if(verbosity>0) cout << "=========== Event " << event->event_number() << " contains " << npart << " particles and " << nvert << " vertices." << endl;

// list of vertices to keep
  vector<HepMC::GenVertex> vkeep;

// trigger
    bool eventok = false;
    for ( HepMC::GenEvent::particle_const_iterator p = event->particles_begin(); p != event->particles_end(); ++p ){
      int pid = (*p)->pdg_id();
      if(abs(pid)==_theTrigger || _theTrigger==0) eventok=true;
    }
    if(!eventok) { return -1; }

// find _theParticle
    for ( HepMC::GenEvent::particle_const_iterator p = event->particles_begin(); p != event->particles_end(); ++p ){
      int pid = (*p)->pdg_id();

      if(abs(pid)==_theParticle) {
           
        bool goodparticle = false;

        // if list of daughters is not empty, check that _theParticle decays to requested daughters
        if(_theDaughters.size()==0) { goodparticle = true; }
        else {
          if ((*p)->end_vertex() ) {
            for ( HepMC::GenVertex::particle_iterator des = (*p)->end_vertex()-> particles_begin(HepMC::descendants);
                                                       des != (*p)->end_vertex()-> particles_end(HepMC::descendants);
                                                     ++des ) {
              for(unsigned int i=0; i<_theDaughters.size(); i++) {
                if(abs((*des)->pdg_id())==_theDaughters[i]) { goodparticle = true; break; }
              }
            }
          }
        }

           if ( goodparticle && (*p)->production_vertex() ) {  // keep this vertex 
                vkeep.push_back(*(*p)->production_vertex());
           }

     } // this is _theParticle
   } // end loop over particles

// loop over vertices and keep only selected ones.
   for ( HepMC::GenEvent::vertex_const_iterator v = event->vertices_begin(); v != event->vertices_end(); ++v ){
     bool goodvertex = false;
     for(unsigned int i=0; i<vkeep.size(); i++) {
       HepMC::GenVertex tmp1 = (*(*v));
       HepMC::GenVertex tmp2 = vkeep[i];
       if(tmp1==tmp2) { goodvertex = true; }
     }
     if(!goodvertex) { event->remove_vertex((*v)); }
   }

// clean up the vertices from particles we don't need
    for ( HepMC::GenEvent::vertex_const_iterator v = event->vertices_begin(); v != event->vertices_end(); ++v ) {
       std::vector<HepMC::GenParticle*> removep;
       for ( HepMC::GenVertex::particle_iterator itpart = (*v)->particles_begin(HepMC::children);
                    itpart != (*v)->particles_end(HepMC::children);
                  ++itpart ) {
         bool keepparticle = false;
         for(unsigned int j=0; j<_theDaughters.size(); j++) { if(abs((*itpart)->pdg_id())==_theDaughters[j] && (*itpart)->status()==1) {keepparticle=true;} }
         if(abs((*itpart)->pdg_id())==_theParticle) {keepparticle=true;}
         if(!keepparticle) { removep.push_back((*itpart)); }
       } // end loop over particles in this vertex
       for(unsigned int k=0; k<removep.size(); k++) { (*v)->remove_particle(removep[k]); }
    } // end look over good vertices

  int partcount=0;
  if(verbosity>0) {
    cout << "FINAL Event " << event->event_number() << " contains " << event->particles_size() << " particles and " << event->vertices_size() << " vertices." << endl;
    cout << "FINAL LIST OF PARTICLES:" << endl;
  }
  for ( HepMC::GenEvent::particle_const_iterator p = event->particles_begin(); p != event->particles_end(); ++p ){
    int pid = (*p)->pdg_id();
    int status = (*p)->status();
    double px = ((*p)->momentum()).px();
    double py = ((*p)->momentum()).py();
    double pz = ((*p)->momentum()).pz();
   double mass  = ((*p)->momentum()).m();
    if(verbosity>0) { cout << pid << " " << mass << " " << status << " " << px << " " << py << " " << pz << endl; }
    partcount++;
  }

  if(partcount==0) { if(verbosity>0) {cout << "EVENT ABORTED: No particles to write out." << endl;} return -1; }
  else { return 0; }

}

void PHHepMCParticleSelectorDecayProduct::SetParticle(const int pid) 
{
  _theParticle = pid;
  return;
}

void PHHepMCParticleSelectorDecayProduct::AddParent(const int pid)
{
  _theParents.push_back(pid);
  return;
}

void PHHepMCParticleSelectorDecayProduct::AddDaughter(const int pid)
{
  _theDaughters.push_back(pid);
  return;
}

