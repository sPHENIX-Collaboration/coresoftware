#ifndef PHGEOMETRY_PHGEOMFILEIMPORT_H
#define PHGEOMETRY_PHGEOMFILEIMPORT_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

//! Fun4All module to import TGeo ROOT Geometry at run time
class PHGeomFileImport : public SubsysReco
{
 public:
  explicit PHGeomFileImport(const std::string &geometry_file);
  ~PHGeomFileImport() override {}

  int InitRun(PHCompositeNode *topNode) override;

 protected:
  std::string m_GeometryFile;
};

#endif
