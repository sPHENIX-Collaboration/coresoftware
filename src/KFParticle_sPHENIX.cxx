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
     initializeBranches();
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
    std::vector<KFParticle> mother, vertex;
    std::vector<std::vector<KFParticle>> daughters, intermediates;
    int nPVs, multiplicity;

    createDecay( topNode, mother, vertex, daughters, intermediates, nPVs, multiplicity );

    if (mother.size() != 0 ) for (unsigned int i = 0; i < mother.size(); ++i) 
    { 
      if ( !m_has_intermediates_sPHENIX ) intermediates.push_back( daughters[i] ); //This is done to avoid a crash, nothing is written to files 
      if ( m_save_output ) fillBranch( topNode, mother[i], vertex[i], daughters[i], intermediates[i], nPVs, multiplicity );
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
