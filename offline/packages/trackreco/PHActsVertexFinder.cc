#include "PHActsVertexFinder.h"

#include <fun4all/Fun4AllReturnCodes.h>





PHActsVertexFinder::PHActsVertexFinder(const std::string &name)
  : PHInitVertexing(name)
{

}


int PHActsVertexFinder::Setup(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsVertexFinder::Process(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}
