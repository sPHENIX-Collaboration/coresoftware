#ifndef KFParticle_DST_H__
#define KFParticle_DST_H__

#include "KFParticle.h"
#include "KFParticle_Container.h"
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                   
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                  
#include <phool/getClass.h>

using namespace std;

class KFParticle;
class KFParticle_Container;
class PHCompositeNode;
class SvtxTrackMap;

class KFParticle_DST
{
  public:

    KFParticle_DST(); //Constructor

    ~KFParticle_DST(); //Destructor

    int createParticleNode(PHCompositeNode* topNode);

    void fillParticleNode(PHCompositeNode* topNode, KFParticle motherParticle,
                          vector<KFParticle> daughters,
                          vector<KFParticle> intermediates);

    void fillParticleNode_Track(PHCompositeNode* topNode, KFParticle motherParticle,
                          vector<KFParticle> daughters,
                          vector<KFParticle> intermediates);

    void fillParticleNode_Particle(PHCompositeNode* topNode, KFParticle motherParticle,
                          vector<KFParticle> daughters,
                          vector<KFParticle> intermediates);

   void printNode(PHCompositeNode* topNode);

protected:

    bool m_has_intermediates_DST;
    bool m_write_track_container;
    bool m_write_particle_container;
    string m_container_name;

private:
    
    SvtxTrackMap* m_recoTrackMap;
    KFParticle_Container* m_recoParticleMap;
    SvtxTrack *m_recoParticle;
    SvtxTrack *buildSvtxTrack( KFParticle particle );

};

#endif //KFParticle_DST_H
