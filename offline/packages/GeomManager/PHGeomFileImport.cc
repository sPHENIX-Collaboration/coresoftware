#include "PHGeomFileImport.h"
#include "PHGeomUtility.h"

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

using namespace std;

PHGeomFileImport::PHGeomFileImport(const std::string & geometry_file) :
    SubsysReco("PHGeomFileImport"), _geometry_file(geometry_file)
{
}

PHGeomFileImport::~PHGeomFileImport()
{
}

int
PHGeomFileImport::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int
PHGeomFileImport::InitRun(PHCompositeNode *topNode)
{
  return PHGeomUtility::ImportGeomFile(topNode, _geometry_file);
}

int
PHGeomFileImport::process_event(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int
PHGeomFileImport::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

