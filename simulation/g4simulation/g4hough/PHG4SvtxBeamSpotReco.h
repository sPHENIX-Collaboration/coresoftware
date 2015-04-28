#ifndef __PHG4SVTXBEAMSPOTRECO__
#define __PHG4SVTXBEAMSPOTRECO__

#include "SvtxVertexMap.h"
#include "SvtxBeamSpot.h"

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

#include <TPrincipal.h>

class PHG4SvtxBeamSpotReco : public SubsysReco {
  
public:

  PHG4SvtxBeamSpotReco(const char * name = "PHG4SvtxBeamSpotReco");
  virtual ~PHG4SvtxBeamSpotReco(){}
  
  //! module initialization
  int Init(PHCompositeNode *topNode){return 0;}
  
  //! run initialization
  int InitRun(PHCompositeNode *topNode);
  
  //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode);
  
private:
  
  // vertex positions
  TPrincipal _pca;

  // pointer to input
  SvtxVertexMap* _vertexes;
  
  // pointer to output
  SvtxBeamSpot* _beamspot;
  
  PHTimeServer::timer _timer;   ///< Timer
};

#endif
