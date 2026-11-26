#ifndef KFPARTICLESPHENIX_KFPARTICLEDST_H
#define KFPARTICLESPHENIX_KFPARTICLEDST_H

#include "KFParticle_Config.h"

#include <KFParticle.h>

#include <string>
#include <vector>

class KFParticle_Container;
class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;

class KFParticle_DST
{
 public:
  /// Constructor. The config is owned by KFParticle_sPHENIX and is not yet constructed when this runs:
  /// store the reference only, never read through it here.
  explicit KFParticle_DST(const KFParticle_Config &config)
    : m_config(config)
  {
  }

  /// Destructor
  virtual ~KFParticle_DST() = default;

  /// Places a KFParticle_Container and SvtxTrackMap on the node tree if they don't exist
  int createParticleNode(PHCompositeNode *topNode);

  /// Simultaneously fills a KFParticle_Container and SvtxTrackMap if they are enabled
  void fillParticleNode(PHCompositeNode *topNode, KFParticle &motherParticle,
                        KFParticle &PV,
                        const std::vector<KFParticle> &daughters,
                        const std::vector<KFParticle> &intermediates);

  /// Called by fillParticleNode, fills an SvtxTrackMap
  void fillParticleNode_Track(PHCompositeNode *topNode, KFParticle &motherParticle,
                              const std::vector<KFParticle> &daughters,
                              const std::vector<KFParticle> &intermediates);

  /// Called by fillParticleNode, fills a KFParitcle_Container
  void fillParticleNode_Particle(PHCompositeNode *topNode, KFParticle &motherParticle,
                                 KFParticle &PV,
                                 std::vector<KFParticle> daughters,
                                 std::vector<KFParticle> intermediates);

  /// Prints contents of KFParticle_Containers and SvtxTrackMaps for an event if they are enabled
  void printNode(PHCompositeNode *topNode);

 private:
  const KFParticle_Config &m_config;
  SvtxTrackMap *m_recoTrackMap = nullptr;
  KFParticle_Container *m_recoParticleMap = nullptr;
  SvtxTrack *buildSvtxTrack(const KFParticle &particle);
};

#endif  // KFPARTICLESPHENIX_KFPARTICLEDST_H
