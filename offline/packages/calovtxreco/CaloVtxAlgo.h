#ifndef CALOVTXALGO_H
#define CALOVTXALGO_H

#include <string>
#include <globalvertex/VertexDefs.h>

class PHCompositeNode;

class CaloVtxAlgo
{
 public:
  virtual ~CaloVtxAlgo() = default;

  // Where set-up occurs:
  //    Calibrations/ML weights etc
  //    Basic information from geometry
  
  virtual int Init(PHCompositeNode * /*topNode*/) { return 0; }

  // Called once per event. Read whatever input this algorithm needs from
  // topNode and fill vtxz [cm]. Return nonzero (e.g.
  // Fun4AllReturnCodes::ABORTEVENT) if no vertex could be found.

  virtual int CalculateVertex(PHCompositeNode *topNode, float &zvtx) = 0;

  // Short identifier used to name this algorithm's output node
  // (e.g. "CALOVTXOUT_MLP") and in log/eval output so results from
  // different algorithms don't collide and can be compared directly.
  virtual std::string Name() const = 0;

  virtual VertexDefs::CALOALGO Algo() const = 0;
};

#endif
