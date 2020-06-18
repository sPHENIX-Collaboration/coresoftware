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

#ifndef KFParticle_sPHENIX_H
#define KFParticle_sPHENIX_H

#include <KFParticle_Tools.h>
#include <KFParticle_nTuple.h>

//KFParticle stuff
#include "KFParticle.h"

//sPHENIX stuff
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/SubsysReco.h>

//ROOT stuff
#include <TFile.h>

class KFParticle_sPHENIX : public SubsysReco, public KFParticle_nTuple, protected KFParticle_Tools
{
 public:
  
  KFParticle_sPHENIX();

  ~KFParticle_sPHENIX();

  int Init(PHCompositeNode *topNode);

  int process_event(PHCompositeNode *topNode); 

  int End(PHCompositeNode *topNode);

  ///Parameters for the user to vary
  void setNumberOfTracks( int num_tracks ) { m_num_tracks = num_tracks; }
  void setMinimumMass( float min_mass ) { m_min_mass = min_mass; }
  void setMaximumMass( float max_mass ) { m_max_mass = max_mass; }
  void setMinimumLifetime( float min_lifetime ) { m_min_lifetime = min_lifetime; }
  void setMaximumLifetime( float max_lifetime ) { m_max_lifetime = max_lifetime; }
  void setMinimumTrackPT( float pt)  { m_track_pt = pt; }
  void setMinimumTrackIPchi2( float ipchi2 ) { m_track_ipchi2 = ipchi2; }
  void setMaximumDaughterDCA( float dca ) { m_comb_DCA = dca; }
  void setFlightDistancechi2( float fdchi2 ) { m_fdchi2 = fdchi2; }
  void setMinDIRA( float dira_min ) { m_dira_min = dira_min; }
  void setMaxDIRA( float dira_max ) { m_dira_max = dira_max; }
  void setMotherPT( float mother_pt ) { m_mother_pt = mother_pt; }
  void setMotherCharge( int mother_charge ) { m_mother_charge = mother_charge; }

  void useMVA( bool require_mva) { m_require_mva = require_mva; }
  void setNumMVAPars( unsigned int nPars ) { m_nPars = nPars; }
  void setMVAVarList( std::string mva_variable_list[ 99 ] ) { for ( int i = 0; i < 99; ++i) m_mva_variable_list[i] = mva_variable_list[i]; }
  void setMVAType( std::string mva_type ) { m_mva_type = mva_type; }
  void setMVAWeightsPath( std::string mva_weights_path ) { m_mva_path = mva_weights_path; }
  void setMVACutValue( float cut_value ) { m_mva_cut_value = cut_value; }

  void setFirstDaughter( std::string name )  { m_daughter_one = name; }
  void setSecondDaughter( std::string name ) { m_daughter_two = name; }
  void setThirdDaughter( std::string name )  { m_daughter_three = name; }
  void setForthDaughter( std::string name ) { m_daughter_four = name; }

  void saveOutput ( bool save ) { m_save_output = save; }
  void setOutputName( std::string name ) { m_outfile_name = name; }

 protected:
 
  bool m_require_mva; 
  bool m_save_output;
  std::string m_outfile_name;
  TFile *m_outfile;

};

#endif //KFParticle_sPHENIX_H

