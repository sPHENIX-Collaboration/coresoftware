/*
 * This file is part of KFParticle package
 * Copyright ( C ) 2007-2019 FIAS Frankfurt Institute for Advanced Studies
 *               2007-2019 Goethe University of Frankfurt
 *               2007-2019 Ivan Kisel <I.Kisel@compeng.uni-frankfurt.de>
 *               2007-2019 Maksym Zyzak
 *
 * KFParticle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * ( at your option ) any later version.
 *
 * KFParticle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/** Mass Hypothesis Codes
 ** These are the codes used in the PDG to identify different particles
 ** that have also been coded into KFParticle
 ** PDG Code | Mass Index | Particle
 **         11 |  0 | Electron
 **         13 |  1 | Muon
 **         19 |  1 | Muon
 **        211 |  2 | Pion ( Charged )
 **        321 |  3 | Kaon ( Charged )
 **       2212 |  4 | Proton
 ** 1000010020 |  5 | Deuteron
 ** 1000010030 |  6 | Triton
 ** 1000020030 |  7 | Helium-3
 ** 1000020040 |  8 | Helium-4
 **       3112 |  9 | Sigma ( - )
 **       3222 | 10 | Sigma ( + )
 **       3312 | 11 | Xi
 **       3334 | 12 | Omega
 ** Use PDG codes in your analysis
 ** All other codes return a pion
 **/

#include "KFParticle_sPHENIX.h"

/// KFParticle constructor
KFParticle_sPHENIX::KFParticle_sPHENIX():
    SubsysReco( "KFPARTICLE" ),
    m_require_mva(false),
    m_save_output(1),
    m_outfile_name("outputData.root")
{}

KFParticle_sPHENIX::~KFParticle_sPHENIX(){} /// KFParticle destructor

int KFParticle_sPHENIX::Init( PHCompositeNode *topNode )
{ 
  if ( m_save_output )
  {
     m_outfile = new TFile(m_outfile_name.c_str(), "RECREATE");
     initializeBranches( m_num_tracks );
  }

  if ( m_require_mva ) 
  {
    TMVA::Reader *reader;
    std::vector<Float_t> MVA_parValues; 
    std::tie( reader, MVA_parValues ) = initMVA();
  }

  return 0;
}

int KFParticle_sPHENIX::process_event( PHCompositeNode *topNode )
{ 
    std::vector<KFParticle> mother, vertex, daughters_1, daughters_2, daughters_3, daughters_4;
    int nPVs, multiplicity;

    createDecay( topNode, mother, vertex, daughters_1, daughters_2, daughters_3, daughters_4, nPVs, multiplicity );

    KFParticle *dummyParticle = new KFParticle();
    for (unsigned int i = 0; i < mother.size(); ++i)
    {
      if (m_num_tracks < 3) daughters_3.push_back(*dummyParticle);
      if (m_num_tracks < 4) daughters_4.push_back(*dummyParticle);
    }


    if (mother.size() != 0 ) for (unsigned int i = 0; i < mother.size(); ++i) 
    { 
      if ( m_save_output ) fillBranch( topNode, mother[i], vertex[i], m_num_tracks, daughters_1[i], daughters_2[i], daughters_3[i], daughters_4[i], nPVs, multiplicity );
    }
    return Fun4AllReturnCodes::EVENT_OK;
}

int KFParticle_sPHENIX::End(PHCompositeNode *topNode)
{
  if ( m_save_output )
  {
     m_outfile->Write();
     m_outfile->Close();
     delete m_outfile;
  }

   return 0;
}
