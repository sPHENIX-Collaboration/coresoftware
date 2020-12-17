#ifndef KFParticle_DST_H__
#define KFParticle_DST_H__

#include "KFParticle_Container.h"
#include "KFParticle_Tools.h"
#include "KFParticle_truthAndDetTools.h"

#include <fun4all/Fun4AllDstInputManager.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrack_v1.h>

#include <KFParticle.h>

using namespace std;

class KFParticle;
class KFParticle_Container;
class PHCompositeNode;
class SvtxTrackMap;

class KFParticle_DST
{
 public:
  ///Constructor
  KFParticle_DST();

  ///Destructor
  ~KFParticle_DST();

  ///Places a KFParticle_Container and SvtxTrackMap on the node tree if they don't exist
  int createParticleNode(PHCompositeNode* topNode);

  ///Simultaneously fills a KFParticle_Container and SvtxTrackMap if they are enabled
  void fillParticleNode(PHCompositeNode* topNode, KFParticle motherParticle,
                        vector<KFParticle> daughters,
                        vector<KFParticle> intermediates);

  ///Called by fillParticleNode, fills an SvtxTrackMap
  void fillParticleNode_Track(PHCompositeNode* topNode, KFParticle motherParticle,
                              vector<KFParticle> daughters,
                              vector<KFParticle> intermediates);

  ///Called by fillParticleNode, fills a KFParitcle_Container
  void fillParticleNode_Particle(PHCompositeNode* topNode, KFParticle motherParticle,
                                 vector<KFParticle> daughters,
                                 vector<KFParticle> intermediates);

  ///Prints contents of KFParticle_Containers and SvtxTrackMaps for an event if they are enabled
  void printNode(PHCompositeNode* topNode);

 protected:
  bool m_has_intermediates_DST;
  bool m_write_track_container;
  bool m_write_particle_container;
  string m_container_name;

 private:
  SvtxTrackMap* m_recoTrackMap = nullptr;
  KFParticle_Container* m_recoParticleMap = nullptr;
  SvtxTrack* m_recoParticle = nullptr;
  SvtxTrack* buildSvtxTrack(KFParticle particle);
};

#endif  //KFParticle_DST_H
