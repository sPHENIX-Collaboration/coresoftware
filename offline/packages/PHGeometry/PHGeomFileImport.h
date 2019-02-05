#ifndef PHGEOMETRY_PHGEOMFILEIMPORT_H
#define PHGEOMETRY_PHGEOMFILEIMPORT_H


#include <fun4all/SubsysReco.h>

#include <ctime>
#include <map>
#include <set>
#include <string>
#include <iostream>

//! Fun4All module to import TGeo ROOT Geometry at run time
class PHGeomFileImport : public SubsysReco
{

public:
  explicit
  PHGeomFileImport(const std::string & geometry_file);
  virtual
  ~PHGeomFileImport();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

protected:

  std::string _geometry_file;

};

#endif
