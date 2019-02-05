#ifndef PHGEOMETRY_PHGEOMFILEIMPORT_H
#define PHGEOMETRY_PHGEOMFILEIMPORT_H

#include <fun4all/SubsysReco.h>

#include <ctime>
#include <iostream>
#include <map>
#include <set>
#include <string>

//! Fun4All module to import TGeo ROOT Geometry at run time
class PHGeomFileImport : public SubsysReco
{
 public:
  explicit PHGeomFileImport(const std::string &geometry_file);
  virtual ~PHGeomFileImport() {}

  int InitRun(PHCompositeNode *topNode);

 protected:
  std::string m_GeometryFile;
};

#endif
