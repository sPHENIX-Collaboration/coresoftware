#ifndef __GLOBALVERTEXRECO_H__
#define __GLOBALVERTEXRECO_H__

//===========================================================
/// \file GlobalVertexReco.h
/// \brief reconstruct the best possible vertexes by
/// combining results from multiple detectors
/// \author Mike McCumber
//===========================================================

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHTimeServer.h>

class PHCompositeNode;

/// \class GlobalVertexReco
///
/// \brief simple truth vertex smearing algorithm
///
class GlobalVertexReco : public SubsysReco {

 public:
 
  GlobalVertexReco(const std::string &name = "GlobalVertexReco");
  virtual ~GlobalVertexReco();
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 private:

  int CreateNodes(PHCompositeNode *topNode);
};

#endif // __GLOBALVERTEXRECO_H__
