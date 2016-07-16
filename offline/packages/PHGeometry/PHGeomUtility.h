#ifndef PHGeomUtility_HH__
#define PHGeomUtility_HH__

#include <ctime>
#include <map>
#include <set>
#include <string>
#include <iostream>

class PHCompositeNode;
class TGeoManager;
class PHGeomTGeo;

//! Toolsets to do geometry operations
class PHGeomUtility
{

protected:

  // static tool sets only
  PHGeomUtility();
  virtual
  ~PHGeomUtility();

public:

  //! DST node name for geometry storage
  static std::string
  GetDSTNodeName()
  {
    return std::string("GEOMETRY");
  }

  //! DST node -> TGeoManager for downstream use
  static TGeoManager *
  GetTGeoManager(PHCompositeNode *topNode);

  //! TGeo ROOT/GDML/Macro file -> DST node with automatic file type discrimination based on file names
  static int
  ImportGeomFile(PHCompositeNode *topNode, const std::string & geometry_file);

//  //! GDML file -> DST node regardless of the file name format
//  static int
//  ImportGDML(PHCompositeNode *topNode, const std::string & gdml_file);

  //! Make a name for tmp geometry file
  //! Geometry files gain a size of ~10MB and it used in translation from Geant4 to DST format.
  //! This tmp file should be on a local file system (/tmp/) and write/deletable
  static std::string GenerateGeometryFileName(const std::string & filename_extension = "gdml");

  //! delete the geometry file after use
  static bool RemoveGeometryFile (const std::string & file_name);


protected:

};

#endif /* PHGeomUtility_HH__ */
