#ifndef FILLTRUTHRECOMATCHMAP_H
#define FILLTRUTHRECOMATCHMAP_H

/**
 * @file trackbase/TrkrMatchDefs.h
 * @author D. Stewart
 * @date February 2023
 * @brief A couple of functions to interpret keys in PHG4ParticleSvtxMap and SvtxPHG4ParticleMap maps keys
 *
 * More detail:
 *  The track matching is filled in EmbRecoMatchContainer objects. For sake of current code practice, these
 *  are interpreted into SvtxPHG4ParticleMap and PHG4ParticleSvtxMap objects. Both those objects are filled with maps
 *  like:
 *    map <track_id{reco,truth}::unsigned int, map< weight::float, track_id{truth,reco}::set<int>>
 *    i.e.   track id of reco or truth, mapped to a second map, which maps weighs to sets of integers on the reco
 *
 *    The weights for a match are given as follows:
 *    Before the decimal, the number of matched clusters;
 *    After the decimal, the number of possible clusters
 *
 *    -------------
 *    Start Example
 *    -------------
 *
 *    An entry in SvtxPHG4ParticleMap could be:
 *      12 -> 41.46 -> { 2, 4 }
 *         -> 18.46 -> { 7 }
 *    which is to say, reco track id 12 has matches to truth tracks 2, 4, and
 *    7. Matches 12->2 and 12->4 weighting key (41.46) indicate that there were
 *    41 matched clusters, that the reco track had 46 clusters. Match 12->7 key
 *    (18.46) indicates that there were 18 matched clusters (and, again, the reco track had 46 clusters)
 *
 *    Assuming that truth tracks 2, 4, and 7, were matched only to reco track 12, and 
 *    each had 45, 44, and 47 clusters, respectively, then the corresponding entries
 *    in PHG4ParticleSvtxMap would be:
 *      2 -> 41.45 { 12 } 
 *      4 -> 41.44 { 12 } 
 *      7 -> 18.47 { 12 }
 *
 */

#include <fun4all/SubsysReco.h> 

class EmbRecoMatchContainer;
class PHCompositeNode;
class PHG4ParticleSvtxMap;
class SvtxPHG4ParticleMap;

class FillTruthRecoMatchMap : public SubsysReco
{
 public:
  FillTruthRecoMatchMap(const std::string &name = "FillTruthRecoMatchMap");

  virtual ~FillTruthRecoMatchMap();

  int Init(PHCompositeNode *) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode * /*topNode*/) override;
  int End(PHCompositeNode *topNode) override;

 private:
  int createNodes(PHCompositeNode *topNode);

  EmbRecoMatchContainer   *m_EmbRecoMatchContainer   {nullptr}; // contianer used to fill the other track matches
  SvtxPHG4ParticleMap     *m_SvtxPHG4ParticleMap     {nullptr}; // reco to truth map, filled for output
  PHG4ParticleSvtxMap     *m_PHG4ParticleSvtxMap     {nullptr}; // truth to reco map, filled for output

};

#endif  // FILLTRUTHRECOMATCHMAP_H
