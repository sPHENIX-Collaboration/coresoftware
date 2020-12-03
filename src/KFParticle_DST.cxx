#include <KFParticle_DST.h>

/*****************/
/* Cameron Dean  */
/*   LANL 2020   */
/* cdean@bnl.gov */ 
/*****************/

/*
 * Class to append reconstructed events to node tree
 * Currently implemented using PHG4Particle
 * A custom class will likely be needed such as PHParticle
 */

//Ideas taken from PHRaveVertexing

KFParticle_DST::KFParticle_DST(){} //Constructor

KFParticle_DST::~KFParticle_DST(){} //Destructor

int KFParticle_DST::createParticleNode(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* lowerNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "Particles"));
  if (!lowerNode)
  {
    lowerNode = new PHCompositeNode("Particles");
    topNode->addNode(lowerNode);
    std::cout<<"Particles node added"<<std::endl;
  }

  m_recoParticleMap = new SvtxTrackMap_v1();
  PHIODataNode<PHObject>* particleNode = new PHIODataNode<PHObject>( m_recoParticleMap, "reconstructedParticles", "PHObject" );
  lowerNode->addNode(particleNode);
  std::cout<<"reconstructedParticles node added"<<std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

void KFParticle_DST::fillParticleNode(PHCompositeNode* topNode, KFParticle motherParticle,
                                      std::vector<KFParticle> daughters,
                                      std::vector<KFParticle> intermediates)
{
  m_recoParticleMap = findNode::getClass<SvtxTrackMap>(topNode, "reconstructedParticles" );  

  //std::string mother_name;
  //if (m_mother_name_DST.empty()) mother_name = "mother";
  //else mother_name = m_mother_name_DST;

  SvtxTrack *m_recoParticle = new SvtxTrack_v1();

  m_recoParticle->set_id(std::abs(motherParticle.GetPDG()));
  m_recoParticle->set_charge(motherParticle.Q());
  m_recoParticle->set_chisq(motherParticle.GetChi2());
  m_recoParticle->set_ndf(motherParticle.GetNDF());
  m_recoParticle->set_x(motherParticle.GetX());
  m_recoParticle->set_y(motherParticle.GetY());
  m_recoParticle->set_z(motherParticle.GetZ());
  m_recoParticle->set_px(motherParticle.GetPx());
  m_recoParticle->set_py(motherParticle.GetPy());
  m_recoParticle->set_pz(motherParticle.GetPz());
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      m_recoParticle->set_error( i, j, motherParticle.GetCovariance( i, j));

  m_recoParticleMap->insert(m_recoParticle);

  std::cout<<"Identifying the particle:"<<std::endl;
  m_recoParticle->identify();
  std::cout<<"Identifying the map:"<<std::endl;
  m_recoParticleMap->identify();

  m_recoParticle->Reset();
 
  if ( m_has_intermediates_DST)
  {
    KFParticle* intermediateArray = &intermediates[0];

    for (unsigned int k = 0; k < intermediates.size(); ++k)
    {
      m_recoParticle->set_id(std::abs(intermediateArray[k].GetPDG()));
      m_recoParticle->set_charge(intermediateArray[k].Q());
      m_recoParticle->set_chisq(intermediateArray[k].GetChi2());
      m_recoParticle->set_ndf(intermediateArray[k].GetNDF());
      m_recoParticle->set_x(intermediateArray[k].GetX());
      m_recoParticle->set_y(intermediateArray[k].GetY());
      m_recoParticle->set_z(intermediateArray[k].GetZ());
      m_recoParticle->set_px(intermediateArray[k].GetPx());
      m_recoParticle->set_py(intermediateArray[k].GetPy());
      m_recoParticle->set_pz(intermediateArray[k].GetPz());
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
          m_recoParticle->set_error( i, j, intermediateArray[k].GetCovariance( i, j));
  
      m_recoParticleMap->insert(m_recoParticle);
      std::cout<<"Identifying the particle:"<<std::endl;
      m_recoParticle->identify();
      std::cout<<"Identifying the map:"<<std::endl;
      m_recoParticleMap->identify();
      m_recoParticle->Reset();
    }
  }

   KFParticle* daughterArray = &daughters[0]; 

  for (unsigned int k = 0; k < daughters.size(); ++k )
  {
    m_recoParticle->set_id(std::abs(daughterArray[k].GetPDG()));
    m_recoParticle->set_charge(daughterArray[k].Q());
    m_recoParticle->set_chisq(daughterArray[k].GetChi2());
    m_recoParticle->set_ndf(daughterArray[k].GetNDF());
    m_recoParticle->set_x(daughterArray[k].GetX());
    m_recoParticle->set_y(daughterArray[k].GetY());
    m_recoParticle->set_z(daughterArray[k].GetZ());
    m_recoParticle->set_px(daughterArray[k].GetPx());
    m_recoParticle->set_py(daughterArray[k].GetPy());
    m_recoParticle->set_pz(daughterArray[k].GetPz());
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 6; ++j)
        m_recoParticle->set_error( i, j, daughterArray[k].GetCovariance( i, j));

    m_recoParticleMap->insert(m_recoParticle);
    std::cout<<"Identifying the particle:"<<std::endl;
    m_recoParticle->identify();
    std::cout<<"Identifying the map:"<<std::endl;
    m_recoParticleMap->identify();
    m_recoParticle->Reset();
  }

}
