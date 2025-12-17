#include "PHGeomFileImport.h"
#include "PHGeomUtility.h"

// PHENIX includes
#include <fun4all/SubsysReco.h>

class PHCompositeNode;

PHGeomFileImport::PHGeomFileImport(const std::string &geometry_file)
  : SubsysReco("PHGeomFileImport")
  , m_GeometryFile(geometry_file)
{
}

int PHGeomFileImport::InitRun(PHCompositeNode *topNode)
{
  return PHGeomUtility::ImportGeomFile(topNode, m_GeometryFile);
}
