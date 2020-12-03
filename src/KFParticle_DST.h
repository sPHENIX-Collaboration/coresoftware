#ifndef KFParticle_DST_H__
#define KFParticle_DST_H__

#include "KFParticle.h"
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

class KFParticle;
class PHCompositeNode;
class SvtxTrackMap;

class KFParticle_DST
{
  public:

    KFParticle_DST(); //Constructor

    ~KFParticle_DST(); //Destructor

    int createParticleNode(PHCompositeNode* topNode);

    void fillParticleNode(PHCompositeNode* topNode, KFParticle motherParticle,
                          std::vector<KFParticle> daughters,
                          std::vector<KFParticle> intermediates);

protected:

    bool m_has_intermediates_DST;

private:
    
    SvtxTrackMap* m_recoParticleMap;
    SvtxTrack *m_recoParticle;
};

#endif //KFParticle_DST_H
