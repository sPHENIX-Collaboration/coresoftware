#include "PHHepMCParticleSelectorB2Jpsi.h"

#include "PHHepMCGenEvent.h"

#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>

using namespace std;

PHHepMCParticleSelectorB2Jpsi::PHHepMCParticleSelectorB2Jpsi(const string &name):
  PHHepMCParticleSelectorBase(name)
{
  _theTrigger = 511;
  _theParticle = 443;
//  _theDaughters.push_back(11);
  _theParents.push_back(445); // also keep B -> chi_c -> J/psi decays
//  _theParents.push_back(511);
//  _theParents.push_back(521);
//  _theParents.push_back(531);
  return;
}

int
PHHepMCParticleSelectorB2Jpsi::InitRun(PHCompositeNode *topNode)
{
  return 0;
}

int
PHHepMCParticleSelectorB2Jpsi::process_event(PHCompositeNode *topNode)
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


// find _thePartcle
    for ( HepMC::GenEvent::particle_const_iterator p = event->particles_begin(); p != event->particles_end(); ++p ){
      int pid = (*p)->pdg_id();
 //     int status = (*p)->status();
 //     double pt = ((*p)->momentum()).perp();
 //     double mass  = ((*p)->momentum()).m();
      if(abs(pid)==_theParticle) {
 //        HepMC::GenVertex* vtxend = (*p)->end_vertex(); // decay vertex
 //        double vx2 = (vtxend->point3d()).x();
 //        double vy2 = (vtxend->point3d()).y();
 //        double vz2 = (vtxend->point3d()).z();

         // production vertex
           bool goodparticle = false;
           if ( (*p)->production_vertex() ) { // just sanity check
             for ( HepMC::GenVertex::particle_iterator mother = (*p)->production_vertex()-> particles_begin(HepMC::parents); 
                                                       mother != (*p)->production_vertex()-> particles_end(HepMC::parents); 
                                                     ++mother ) {
                //HepMC::GenVertex* tmpv0 = (*p)->production_vertex();
                //cout << "      J/psi production vertex: " << tmpv0->point3d().x() << " " << tmpv0->point3d().y() << " " << tmpv0->point3d().z() << endl; 
                //std::cout << "      parent: " << (*mother)->pdg_id() << endl;
                for(unsigned int i=0; i<_theParents.size(); i++) {
                  if(abs((*mother)->pdg_id())==_theParents[i]) { 
                    //HepMC::GenVertex* tmpv = (*mother)->production_vertex();
                    //cout << "         parent production vertex: " << tmpv->point3d().x() << " " << tmpv->point3d().y() << " " << tmpv->point3d().z() << endl; 
                    vkeep.push_back(*(*p)->production_vertex());
                    vkeep.push_back(*(*mother)->production_vertex());
                    goodparticle = true;
                    break;
                  }
                }
             }
           }

         // decay vertex
           if ( goodparticle && (*p)->end_vertex() ) {
             for ( HepMC::GenVertex::particle_iterator des = (*p)->end_vertex()-> particles_begin(HepMC::descendants);
                                                       des != (*p)->end_vertex()-> particles_end(HepMC::descendants);
                                                     ++des ) {
                //HepMC::GenVertex* tmpv0 = (*p)->end_vertex();
                //cout << "      J/psi decay vertex: " << tmpv0->point3d().x() << " " << tmpv0->point3d().y() << " " << tmpv0->point3d().z() << endl; 
                //std::cout << "      descendant: " << (*des)->pdg_id() << endl;
                for(unsigned int i=0; i<_theDaughters.size(); i++) {
                  if(abs((*des)->pdg_id())==_theDaughters[i]) { 
                    //HepMC::GenVertex* tmpv = (*des)->production_vertex();
                    //cout << "         descendant production vertex: " << tmpv->point3d().x() << " " << tmpv->point3d().y() << " " << tmpv->point3d().z() << endl; 
                    vkeep.push_back(*(*p)->end_vertex());
                    break;
                  }
                }
             }
           }

      } // this is _theParticle
    } // end loop over particles
    //cout << "         Found " << vkeep.size() << " vertices to keep." << endl;

// loop over vertices and keep only selected ones.
    //cout << "looping over vertices..." << endl;
    for ( HepMC::GenEvent::vertex_const_iterator v = event->vertices_begin(); v != event->vertices_end(); ++v ){
      bool goodvertex = false;
      for(unsigned int i=0; i<vkeep.size(); i++) {
        HepMC::GenVertex tmp1 = (*(*v));
        HepMC::GenVertex tmp2 = vkeep[i];
        if(tmp1==tmp2) { goodvertex = true; } 
      } 
      if(!goodvertex) { event->remove_vertex((*v)); }
    }    

/*
if(verbosity) {
cout << "INTERMEDIATE Event " << event->event_number() << " contains " << event->particles_size() << " particles and " << event->vertices_size() << " vertices." << endl;
cout << "INTERMEDIATE LIST OF PARTICLES:" << endl;
    for ( HepMC::GenEvent::particle_const_iterator p = event->particles_begin(); p != event->particles_end(); ++p ){
      int pid = (*p)->pdg_id();
      int status = (*p)->status();
      double pt = ((*p)->momentum()).perp();
      double mass  = ((*p)->momentum()).m();
      HepMC::GenVertex* tmpv = (*p)->production_vertex();
      cout << pid << " " << mass << " " << status << " " << pt << " " << tmpv << endl;
    }
}
*/

// clean up the vertices
//cout << "cleaning up the vertices..." << endl;
    for ( HepMC::GenEvent::vertex_const_iterator v = event->vertices_begin(); v != event->vertices_end(); ++v ) {
       //cout << "vertex # " << vtxcount << "        " << (*v)->point3d().x() << " " << (*v)->point3d().y() << " " << (*v)->point3d().z() << endl;
       std::vector<HepMC::GenParticle*> removep;
       for ( HepMC::GenVertex::particle_iterator itpart = (*v)->particles_begin(HepMC::children);
                    itpart != (*v)->particles_end(HepMC::children);
                  ++itpart ) {
         bool keepparticle = false;
         for(unsigned int j=0; j<_theDaughters.size(); j++) { if(abs((*itpart)->pdg_id())==_theDaughters[j] && (*itpart)->status()==1) {keepparticle=true;} }
         for(unsigned int j=0; j<_theParents.size(); j++) { if(abs((*itpart)->pdg_id())==_theParents[j] && (*itpart)->status()==1) {keepparticle=true;} }
         if(abs((*itpart)->pdg_id())==_theParticle) {keepparticle=true;}
         if(!keepparticle) { removep.push_back((*itpart)); }
       } // end loop over particles in this vertex
       //for(unsigned int k=0; k<removep.size(); k++) { HepMC::GenParticle *tmp = (*v)->remove_particle(removep[k]); }
       for(unsigned int k=0; k<removep.size(); k++) { (*v)->remove_particle(removep[k]); }
    }

/*
    bool abortevent = false;
    for ( HepMC::GenEvent::particle_const_iterator p = event->particles_begin(); p != event->particles_end(); ++p ){
      int status = (*p)->status();
      if(status!=1) abortevent = true;
    }
*/

int partcount=0;
//if(verbosity && !abortevent) {
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
    //double pt = ((*p)->momentum()).perp();
    double mass  = ((*p)->momentum()).m();
    if(verbosity>0) { cout << pid << " " << mass << " " << status << " " << px << " " << py << " " << pz << endl; }
    partcount++;
  }

  if(partcount==0) { if(verbosity>0) {cout << "EVENT ABORTED: No particles to write out." << endl;} return -1; }
  else { return 0; }
//  if(abortevent) { cout << "EVENT ABORTED." << endl; return -1; }
//  else { return 0; }
//  return 0;

}



