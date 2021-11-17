/**
 * This file is part of KFParticle package
 * Copyright (C) 2007-2019 FIAS Frankfurt Institute for Advanced Studies
 *               2007-2019 Goethe University of Frankfurt
 *               2007-2019 Ivan Kisel <I.Kisel@compeng.uni-frankfurt.de>
 *               2007-2019 Maksym Zyzak
 *
 * KFParticle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KFParticle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef KFPARTICLESPHENIX_KFPARTICLEPARTICLELIST_H
#define KFPARTICLESPHENIX_KFPARTICLEPARTICLELIST_H

#include <map>

typedef std::pair<int, float> particle_pair;

/**
 *  This class contains the information for each type of particle
 *  A map is created with the particle name as the key. 
 *  This will then return the PDG ID and mass in GeV.
 *  The map element format is, for example:
 *  particleMasses["Bs0"] = make_pair(531, 5.366);
 */
class KFParticleDatabase;

class KFParticle_particleList
{
 public:
  KFParticle_particleList(){}

  virtual ~KFParticle_particleList(){}

  std::map<std::string, particle_pair> getParticleList();

  float returnPDGMass(const int pdgIndex);
};

#endif  //KFPARTICLESPHENIX_KFPARTICLEPARTICLELIST_H
