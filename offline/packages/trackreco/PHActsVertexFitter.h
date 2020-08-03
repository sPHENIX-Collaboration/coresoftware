#ifndef TRACKRECO_PHACTSVERTEXFITTER_H
#define TRACKRECO_PHACTSVERTEXFITTER_H

#include <fun4all/SubsysReco.h>
#include "PHActsSourceLinks.h"

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;



class PHActsVertexFitter : public SubsysReco
{
 public:
  PHActsVertexFitter();
  virtual ~PHActsVertexFitter(){}
  int process_event(PHCompositeNode *topNode);
  int Init(PHCompositeNode *topNode);
  int End (PHCompositeNode *topNode);


 private:

  int getNodes(PHCompositeNode *topNode);
  
  SvtxTrackMap *m_svtxTrackMap;
  ActsTrackingGeometry *m_tGeometry;
};

#endif //TRACKRECO_PHACTSVERTEXFITTER_H 
