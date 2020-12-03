/*
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

#ifndef KFParticle_particleList_H
#define KFParticle_particleList_H

#include <map>

using namespace std;

typedef pair<int, float> particle_pair;

class KFParticleDatabase;

class KFParticle_particleList  
{
 public:
  
  KFParticle_particleList();

  ~KFParticle_particleList();

  //map<string, float> getParticleList();
  map<string, particle_pair> getParticleList();

  float returnPDGMass( const int pdgIndex);
};

#endif //KFParticle_particleList_H
