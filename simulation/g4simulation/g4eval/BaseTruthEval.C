#include "BaseTruthEval.h"

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <cassert>
#include <iostream>

using namespace std;

BaseTruthEval::BaseTruthEval(PHCompositeNode* topNode)
  : _truthinfo(nullptr),
    _strict(false),
    _verbosity(1),
    _errors(0) {
  get_node_pointers(topNode);
}

BaseTruthEval::~BaseTruthEval() {
  if (_verbosity > 0) {
    if ((_errors > 0)||(_verbosity > 1)) {
      cout << "BaseTruthEval::~BaseTruthEval() - Error Count: " << _errors << endl;
    }
  }
}

void BaseTruthEval::next_event(PHCompositeNode* topNode) {
  get_node_pointers(topNode);
}

int BaseTruthEval::get_embed(PHG4Particle* particle) {

  if (!has_reduced_node_pointers()) {++_errors; return 0;}
  
  if (_strict) {assert(particle);}
  else if (!particle) {++_errors; return 0;}

//  if (!is_primary(particle)) return 0;
  // allow the secondary particles to be tagged with its embedding ID of tis primary particle

  PHG4Particle* primary = get_primary_particle(particle);
  if (_strict) {assert(primary);}
  else if (!primary) {++_errors; return 0;}
  
  return _truthinfo->isEmbeded(primary->get_track_id());
}

PHG4VtxPoint* BaseTruthEval::get_vertex(PHG4Particle* particle) {

  if (!has_reduced_node_pointers()) {++_errors; return nullptr;} 
  
  if (_strict) {assert(particle);}
  else if (!particle) {++_errors; return nullptr;}

  PHG4VtxPoint* vtx = _truthinfo->GetVtx( particle->get_vtx_id() );
  if (_strict) {assert(vtx);}
  else if (!vtx) {_errors++;}
 
  return vtx;
}

bool BaseTruthEval::is_primary(PHG4Shower* shower) {

  if (!has_reduced_node_pointers()) {++_errors; return false;}
  
  if (_strict) {assert(shower);}
  else if (!shower) {++_errors; return false;}
  
  bool is_primary = false;
  if (shower->get_parent_shower_id() == 0) {
    is_primary = true;
  }
  
  return is_primary;
}

bool BaseTruthEval::is_primary(PHG4Particle* particle) {

  if (!has_reduced_node_pointers()) {++_errors; return false;}
  
  if (_strict) {assert(particle);}
  else if (!particle) {++_errors; return false;}
  
  bool is_primary = false;
  if (particle->get_parent_id() == 0) {
    is_primary = true;
  }
  
  return is_primary;
}

PHG4Shower* BaseTruthEval::get_primary_shower(PHG4Shower* shower) {

  if (!has_reduced_node_pointers()) {++_errors; return nullptr;}
  
  if (_strict) {assert(shower);}
  else if (!shower) {++_errors; return nullptr;}

  if (is_primary(shower)) return shower;

  while (!is_primary(shower)) {
    shower = _truthinfo->GetShower( shower->get_parent_shower_id() );

    if (_strict) {assert(shower);}
    else if (!shower) {++_errors; break;}
  }

  return shower;
}

PHG4Shower* BaseTruthEval::get_primary_shower(PHG4Particle* particle) {

  if (!has_reduced_node_pointers()) {++_errors; return nullptr;}
  
  if (_strict) {assert(particle);}
  else if (!particle) {++_errors; return nullptr;}

  if (!is_primary(particle)) particle = get_primary_particle(particle);
  
  PHG4Shower* returnval = nullptr;
  
  PHG4TruthInfoContainer::ShowerRange range = _truthinfo->GetPrimaryShowerRange();
  for (PHG4TruthInfoContainer::ShowerIterator iter = range.first;
       iter != range.second;
       ++iter) {
    PHG4Shower* shower = iter->second;
    if (shower->get_parent_particle_id() == particle->get_track_id()) {
      returnval = shower;
      break;
    }
  }

  return returnval;
}

PHG4Particle* BaseTruthEval::get_primary_particle(PHG4Particle* particle) {

  if (!has_reduced_node_pointers()) {++_errors; return nullptr;}
  
  if (_strict) {assert(particle);}
  else if (!particle) {++_errors; return nullptr;}

  if (is_primary(particle)) return particle;
  
  PHG4Particle* returnval = _truthinfo->GetPrimaryParticle( particle->get_primary_id() );

  if (_strict) {assert(returnval);}
  else if (!returnval) {++_errors;}
  
  return returnval;
}

PHG4Particle* BaseTruthEval::get_primary_particle(PHG4Shower* shower) {

  if (!has_reduced_node_pointers()) {++_errors; return nullptr;}
  
  if (_strict) {assert(shower);}
  else if (!shower) {++_errors; return nullptr;}

  PHG4Particle* parent_particle = _truthinfo->GetParticle( shower->get_parent_particle_id() );

  if (_strict) {assert(parent_particle);}
  else if (!parent_particle) {++_errors; return nullptr;}
  
  PHG4Particle* primary_particle = get_primary_particle(parent_particle);

  if (_strict) {assert(primary_particle);}
  else if (!primary_particle) {++_errors; return nullptr;}
    
  return primary_particle;
}

std::set<PHG4Shower*> BaseTruthEval::all_secondary_showers(PHG4Shower* shower) {

  if (!has_reduced_node_pointers()) {++_errors; return std::set<PHG4Shower*>();}
  
  if (_strict) {assert(shower);}
  else if (!shower) {++_errors; return std::set<PHG4Shower*>();}
  
  std::set<PHG4Shower*> subshowers;
  
  PHG4TruthInfoContainer::ShowerRange range = _truthinfo->GetSecondaryShowerRange();
  for (PHG4TruthInfoContainer::ShowerIterator iter = range.first;
       iter != range.second;
       ++iter) {
    PHG4Shower* shower = iter->second;

    if (_strict) {assert(shower);}
    else if (!shower) {++_errors;}

    if (shower) {
      if (shower->get_parent_shower_id() == shower->get_id()) {
	subshowers.insert(shower);
      }
    }
  }
  
  return subshowers;
}

bool BaseTruthEval::are_same_shower(PHG4Shower* s1, PHG4Shower* s2) {

  if (!has_reduced_node_pointers()) {++_errors; return false;}
  
  if (_strict) {
    assert(s1);
    assert(s2);
  } else if (!s1||!s2) {
    ++_errors;
    return false;
  }

  if (s1->get_id() == s2->get_id()) return true;
  return false;
}

bool BaseTruthEval::are_same_particle(PHG4Particle* p1, PHG4Particle* p2) {

  if (!has_reduced_node_pointers()) {++_errors; return false;}
  
  if (_strict) {
    assert(p1);
    assert(p2);
  } else if (!p1||!p2) {
    ++_errors;
    return false;
  }

  if (p1->get_track_id() == p2->get_track_id()) return true;
  return false;
}

bool BaseTruthEval::are_same_vertex(PHG4VtxPoint* vtx1, PHG4VtxPoint* vtx2) {

  if (!has_reduced_node_pointers()) {++_errors; return false;}
  
  if (_strict) {
    assert(vtx1);
    assert(vtx2);
  } else if (!vtx1||!vtx2) {
    ++_errors;
    return false;
  }

  if (vtx1->get_id() == vtx2->get_id()) return true;
  return false;
}

PHG4Particle* BaseTruthEval::get_particle(PHG4Hit* g4hit) {

  if (!has_reduced_node_pointers()) {++_errors; return nullptr;}
  
  if (_strict) {assert(g4hit);}
  else if (!g4hit) {++_errors; return nullptr;}
  
  PHG4Particle* particle = _truthinfo->GetParticle( g4hit->get_trkid() );
  if (_strict) {assert(particle);}
  else if (!particle) {++_errors;}
  
  return particle;
}

PHG4Shower* BaseTruthEval::get_primary_shower(PHG4Hit* g4hit) {

  if (!has_reduced_node_pointers()) {++_errors; return nullptr;}
  
  if (_strict) {assert(g4hit);}
  else if (!g4hit) {++_errors; return nullptr;}
  
  PHG4Shower* shower = _truthinfo->GetShower( g4hit->get_shower_id() );
  if (_strict) {assert(shower);}
  else if (!shower) {++_errors;}
  
  return shower;
}

PHG4Particle* BaseTruthEval::get_primary_particle(PHG4Hit* g4hit) {

  if (!has_reduced_node_pointers()) {++_errors; return nullptr;}
  
  if (_strict) {assert(g4hit);}
  else if (!g4hit) {++_errors; return nullptr;}

  PHG4Particle* particle = get_particle(g4hit);
  PHG4Particle* primary = get_primary_particle(particle);

  if (_strict) {assert(primary);}
  else if (!primary) {++_errors;}
  
  return primary;
}

bool BaseTruthEval::is_g4hit_from_primary_shower(PHG4Hit* g4hit, PHG4Shower* shower) {

  if (!has_reduced_node_pointers()) {++_errors; return false;}
   
  if (_strict) {
    assert(g4hit);
    assert(shower);
  } else if (!g4hit||!shower) {
    ++_errors;
    return false;
  }

  if (g4hit->get_shower_id() == shower->get_id()) {
    return true;    
  }

  return false;
}

bool BaseTruthEval::is_g4hit_from_particle(PHG4Hit* g4hit, PHG4Particle* particle) {

  if (!has_reduced_node_pointers()) {++_errors; return false;}
   
  if (_strict) {
    assert(g4hit);
    assert(particle);
  } else if (!g4hit||!particle) {
    ++_errors;
    return false;
  }

  if (g4hit->get_trkid() == particle->get_track_id()) {
    return true;    
  }

  return false;
}

void BaseTruthEval::get_node_pointers(PHCompositeNode* topNode) {

  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!_truthinfo) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }
  
  return;
}

bool BaseTruthEval::has_reduced_node_pointers() {

  if (_strict) assert(_truthinfo);
  else if (!_truthinfo) return false;
  
  return true;
}
