#include <KFParticle_DST.h>
#include <KFParticle_Tools.h>
#include <KFParticle_truthAndDetTools.h>

/*****************/
/* Cameron Dean  */
/*   LANL 2020   */
/* cdean@bnl.gov */ 
/*****************/

/*
 * Class to append reconstructed events to node tree
 */

//Ideas taken from PHRaveVertexing

using namespace std;

KFParticle_Tools kfpTupleTools_DST;
KFParticle_truthAndDetTools kfpTruthTools_DST;

KFParticle_DST::KFParticle_DST():
  m_write_track_container(true),
  m_write_particle_container(true)
{} //Constructor

KFParticle_DST::~KFParticle_DST(){} //Destructor

int KFParticle_DST::createParticleNode(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* lowerNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  //PHCompositeNode* lowerNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "Particles"));
  if (!lowerNode)
  {
    lowerNode = new PHCompositeNode("DST");
    //lowerNode = new PHCompositeNode("Particles");
    topNode->addNode(lowerNode);
    cout<<"Particles node added"<<endl;
  }

  string baseName, trackNodeName, particleNodeName;
  if (m_container_name.empty()) baseName = "reconstructedParticles";
  else baseName = m_container_name;
  //Cant have forward slashes in DST or else you make a subdirectory on save!!!
  string fwd_slsh = "/", undrscr = "_"; size_t pos;
  while ((pos = baseName.find(fwd_slsh)) != string::npos) baseName.replace(pos, 1, undrscr);

  trackNodeName = baseName + "_SvtxTrackMap";
  particleNodeName = baseName + "_KFParticle_Container";

  if ( m_write_track_container )
  {
    m_recoTrackMap = new SvtxTrackMap_v1();
    PHIODataNode<PHObject>* trackNode = new PHIODataNode<PHObject>( m_recoTrackMap, trackNodeName.c_str(), "PHObject" );
    lowerNode->addNode(trackNode);
    printf("%s node added\n", trackNodeName.c_str());
  }

  if ( m_write_particle_container )
  {
    m_recoParticleMap = new KFParticle_Container();
    PHIODataNode<PHObject>* particleNode = new PHIODataNode<PHObject>( m_recoParticleMap, particleNodeName.c_str(), "PHObject" );
    lowerNode->addNode(particleNode);
    printf("%s node added\n", particleNodeName.c_str());
  }

  if ( !m_write_track_container && !m_write_particle_container )
  {
    cout<<"You have asked to put your selection on the node tree but disabled both the SvtxTrackMap and KFParticle_Container\n";
    cout<<"Check your options"<<endl;
    exit(0);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


void KFParticle_DST::fillParticleNode(PHCompositeNode* topNode, KFParticle motherParticle,
                                      vector<KFParticle> daughters,
                                      vector<KFParticle> intermediates)
{
  if ( m_write_track_container ) fillParticleNode_Track(topNode, motherParticle,
                                                        daughters, intermediates);

  if ( m_write_particle_container ) fillParticleNode_Particle(topNode, motherParticle,
                                                              daughters, intermediates);
}

void KFParticle_DST::fillParticleNode_Track(PHCompositeNode* topNode, KFParticle motherParticle,
                                      vector<KFParticle> daughters,
                                      vector<KFParticle> intermediates)
{
  string baseName, trackNodeName;
  if (m_container_name.empty()) baseName = "reconstructedParticles";
  else baseName = m_container_name;
  //Cant have forward slashes in DST or else you make a subdirectory on save!!!
  string fwd_slsh = "/", undrscr = "_"; size_t pos;
  while ((pos = baseName.find(fwd_slsh)) != string::npos) baseName.replace(pos, 1, undrscr);

  trackNodeName = baseName + "_SvtxTrackMap";

  m_recoTrackMap = findNode::getClass<SvtxTrackMap>(topNode, trackNodeName.c_str() );  

  SvtxTrack *m_recoTrack = new SvtxTrack_v1();

  m_recoTrack = buildSvtxTrack(motherParticle);
  m_recoTrackMap->insert(m_recoTrack);
  m_recoTrack->Reset();
 
  if ( m_has_intermediates_DST)
  {
    KFParticle* intermediateArray = &intermediates[0];

    for (unsigned int k = 0; k < intermediates.size(); ++k)
    {
      m_recoTrack = buildSvtxTrack(intermediateArray[k]);
      m_recoTrackMap->insert(m_recoTrack);
      m_recoTrack->Reset();
    }
  }

  SvtxTrackMap* originalTrackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap" );
  SvtxTrackMap* originalTrackMap_copy = new SvtxTrackMap(*originalTrackMap); //Copy needed to avoid memory corruption when using getTrack here and in truth matching
  KFParticle* daughterArray = &daughters[0]; 

  for (unsigned int k = 0; k < daughters.size(); ++k )
  {

    if ( originalTrackMap_copy->size() == 0 ) 
    { 
      cout << "There was no orginal track map found, the tracks will have no cluster information!" << endl;   
      m_recoTrack = buildSvtxTrack(daughterArray[k]);
    }
    else m_recoTrack = kfpTruthTools_DST.getTrack(daughterArray[k].Id(), originalTrackMap_copy);

    m_recoTrackMap->insert(m_recoTrack);
    m_recoTrack->Reset();
  }

}

void KFParticle_DST::fillParticleNode_Particle(PHCompositeNode* topNode, KFParticle motherParticle,
                                      vector<KFParticle> daughters,
                                      vector<KFParticle> intermediates)
{   
  string baseName, particleNodeName;
  if (m_container_name.empty()) baseName = "reconstructedParticles";
  else baseName = m_container_name;
  //Cant have forward slashes in DST or else you make a subdirectory on save!!!
  string fwd_slsh = "/", undrscr = "_"; size_t pos;
  while ((pos = baseName.find(fwd_slsh)) != string::npos) baseName.replace(pos, 1, undrscr);

  particleNodeName = baseName + "_KFParticle_Container";

  m_recoParticleMap = findNode::getClass<KFParticle_Container>(topNode, particleNodeName.c_str() );  

  m_recoParticleMap->insert(&motherParticle);

  if ( m_has_intermediates_DST)
  {
    KFParticle* intermediateArray = &intermediates[0];
    for (unsigned int k = 0; k < intermediates.size(); ++k)
      m_recoParticleMap->insert(&intermediateArray[k]);
  }

  KFParticle* daughterArray = &daughters[0]; 
  for (unsigned int k = 0; k < daughters.size(); ++k )
      m_recoParticleMap->insert(&daughterArray[k]);
}

SvtxTrack* KFParticle_DST::buildSvtxTrack( KFParticle particle )
{
  SvtxTrack *track = new SvtxTrack_v1();

  track->set_id(abs(particle.GetPDG()));
  track->set_charge((int) particle.GetQ());
  track->set_chisq(particle.GetChi2());
  track->set_ndf(particle.GetNDF());

  track->set_x(particle.GetX());
  track->set_y(particle.GetY());
  track->set_z(particle.GetZ());

  track->set_px(particle.GetPx());
  track->set_py(particle.GetPy());
  track->set_pz(particle.GetPz());

  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      track->set_error( i, j, particle.GetCovariance( i, j));

  return track;
}

void KFParticle_DST::printNode(PHCompositeNode* topNode)
{
  string baseName, trackNodeName, particleNodeName;
  if (m_container_name.empty()) baseName = "reconstructedParticles";
  else baseName = m_container_name;
  //Cant have forward slashes in DST or else you make a subdirectory on save!!!
  string fwd_slsh = "/", undrscr = "_"; size_t pos;
  while ((pos = baseName.find(fwd_slsh)) != string::npos) baseName.replace(pos, 1, undrscr);


  if ( m_write_track_container )
  {
    trackNodeName = baseName + "_SvtxTrackMap";
    SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, trackNodeName.c_str());  
    for ( SvtxTrackMap::Iter iter = trackmap->begin(); iter != trackmap->end(); ++iter )
    {  
       SvtxTrack *track = iter->second;
       track->identify();
    }
  }

  if ( m_write_particle_container )
  {
    particleNodeName = baseName + "_KFParticle_Container";
    KFParticle_Container *particlemap = findNode::getClass<KFParticle_Container>(topNode, particleNodeName.c_str());  
    for ( KFParticle_Container::Iter iter = particlemap->begin(); iter != particlemap->end(); ++iter )
    {
      KFParticle *particle = iter->second;
      kfpTupleTools_DST.identify(*particle);
    }
  }
}
