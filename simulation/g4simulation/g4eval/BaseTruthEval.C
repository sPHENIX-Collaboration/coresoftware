
#include "BaseTruthEval.h"

#include <fun4all/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Hit.h>

#include <cstdlib>
#include <float.h>
#include <cassert>
#include <iostream>

using namespace std;

BaseTruthEval::BaseTruthEval(PHCompositeNode* topNode)
  : _truthinfo(NULL),
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
      
PHG4Particle* BaseTruthEval::get_particle(PHG4Hit* g4hit) {

  if (_strict) {assert(g4hit);}
  else if (!g4hit) {++_errors; return NULL;}
  
  PHG4Particle* particle = _truthinfo->GetHit( g4hit->get_trkid() );
  if (_strict) {assert(particle);}
  else if (!particle) {++_errors;}
  
  return particle;
}

int BaseTruthEval::get_embed(PHG4Particle* particle) {

  if (_strict) {assert(particle);}
  else if (!particle) {++_errors; return 0;}

  if (!is_primary(particle)) return 0;

  PHG4Particle* primary = get_primary(particle);
  if (_strict) {assert(primary);}
  else if (!primary) {++_errors; return 0;}
  
  return _truthinfo->isEmbeded(primary->get_track_id());
}

PHG4VtxPoint* BaseTruthEval::get_vertex(PHG4Particle* particle) {

  if (_strict) {assert(particle);}
  else if (!particle) {++_errors; return NULL;}

  if (particle->get_primary_id() == -1) {
    return _truthinfo->GetPrimaryVtx( particle->get_vtx_id() );  
  }

  PHG4VtxPoint* vtx = _truthinfo->GetVtx( particle->get_vtx_id() );
  if (_strict) {assert(vtx);}
  else if (!vtx) {_errors++;}
 
  return vtx;
}

bool BaseTruthEval::is_primary(PHG4Particle* particle) {

  if (_strict) {assert(particle);}
  else if (!particle) {return false;}
  
  bool is_primary = false;
  if (particle->get_parent_id() == 0) {
    is_primary = true;
  }
  
  return is_primary;
}

PHG4Particle* BaseTruthEval::get_primary(PHG4Hit* g4hit) {

  if (_strict) {assert(g4hit);}
  else if (!g4hit) {++_errors; return NULL;}

  PHG4Particle* particle = get_particle(g4hit);
  PHG4Particle* primary = get_primary(particle);

  if (_strict) {assert(primary);}
  else if (!primary) {++_errors;}
  
  return primary;
}

PHG4Particle* BaseTruthEval::get_primary(PHG4Particle* particle) {

  if (_strict) {assert(particle);}
  else if (!particle) {++_errors; return NULL;}

  // always report the primary from the Primary Map regardless if a
  // primary from the full Map was the argument
  PHG4Particle* returnval = NULL;
  if (particle->get_primary_id() != -1) {
    returnval = _truthinfo->GetPrimaryHit( particle->get_primary_id() );
  } else {
    returnval = _truthinfo->GetPrimaryHit( particle->get_track_id() );
  }

  if (_strict) {assert(returnval);}
  else if (!returnval) {++_errors;}
  
  return returnval;
}

 bool BaseTruthEval::is_g4hit_from_particle(PHG4Hit* g4hit, PHG4Particle* particle) {

  if (_strict) {
    assert(g4hit);
    assert(particle);
  } else if (!g4hit||!particle) {
    ++_errors;
    return false;
  }

  bool is_from_particle = false;

  PHG4Particle* candidate = get_particle(g4hit);
  if (_strict) {assert(candidate);}
  else if (!candidate) {++_errors;}

  if (candidate) {
    if (are_same_particle(candidate,particle)) {
      is_from_particle = true;
    }
  }

  return is_from_particle;
}

bool BaseTruthEval::are_same_particle(PHG4Particle* p1, PHG4Particle* p2) {

  if (_strict) {
    assert(p1);
    assert(p2);
  } else if (!p1||!p2) {
    ++_errors;
    return false;
  }

  bool are_same_particle = false;

  // get ready for some insanity... I have to deal with the case where
  // the primary particles are copied between the two truth containers
  // but have different track ids
  
  if ((p1->get_primary_id() == -1) &&
      (p2->get_primary_id() == -1)) {
    // both particles are from the primary truth map, use track ids
    if (p1->get_track_id() == p2->get_track_id()) {
      are_same_particle = true;
    }
  } else if ((p1->get_primary_id() != -1) &&
	     (p2->get_primary_id() != -1)) {
    // both particles are from the total truth map, use track ids
    if (p1->get_track_id() == p2->get_track_id()) {
      are_same_particle = true;
    }
  } else if (is_primary(p1) && is_primary(p2)) {
    // both particles are primary
    if ((p1->get_primary_id() == -1) &&
	(p2->get_primary_id() != -1)) {
      // only the first particle is from the primary truth map, use one track id and one primary id
      if (p1->get_track_id() == p2->get_primary_id()) {
	are_same_particle = true;
      }
    } else if ((p1->get_primary_id() != -1) &&
	       (p2->get_primary_id() == -1)) {
      // only the second particle is from the primary truth map, use one track id and one primary id
      if (p1->get_primary_id() == p2->get_track_id()) {
	are_same_particle = true;
      }
    }
  }
 
  return are_same_particle;
}

bool BaseTruthEval::are_same_vertex(PHG4VtxPoint* vtx1, PHG4VtxPoint* vtx2) {

  if (_strict) {
    assert(vtx1);
    assert(vtx2);
  } else if (!vtx1||!vtx2) {
    ++_errors;
    return false;
  }

  bool are_same_vertex = false;

  // I will need to deal with copies of vertex points too, argh!
  // Here I have very few options, we can either match the vertexes by
  // id or by position there is no ancestry to utilize

  // in the embedded dst it is 100% sure the first method will fail
  // since I can't trace between vertexes and the vertex ids have changed
  // I guess I will have to resort to the second method, which means there
  // is some possiblity of having a position collision between two vertexes
  // that are in fact different in terms of content
  
  //  if (vtx1->get_id() == vtx2->get_id()) {
  //    are_same_vertex = true;
  //  }

  if ((vtx1->get_x() == vtx2->get_x()) &&
      (vtx1->get_y() == vtx2->get_y()) &&
      (vtx1->get_z() == vtx2->get_z()) &&
      (vtx1->get_t() == vtx2->get_t())) {
    are_same_vertex = true;
  }
 
  return are_same_vertex;
}
  
void BaseTruthEval::get_node_pointers(PHCompositeNode* topNode) {

  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!_truthinfo) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }
  
  return;
}
