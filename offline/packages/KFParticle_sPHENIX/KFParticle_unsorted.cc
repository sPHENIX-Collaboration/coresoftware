// TEMPORARY: see KFParticle_unsorted.h. The function bodies are moved from KFParticle_Tools.cc unchanged.

#include "KFParticle_unsorted.h"

#include <TDatabasePDG.h>

#include <iostream>  // for operator<<, basic_ostream

bool KFParticle_unsorted::findParticle(const std::string &particle)
{
  bool particleFound = true;
  if (!TDatabasePDG::Instance()->GetParticle(particle.c_str()))
  {
    std::cout << "The particle, " << particle << " is not recognized" << std::endl;
    std::cout << "Check TDatabasePDG for a list of available particles" << std::endl;
    particleFound = false;
  }

  return particleFound;
}

