#ifndef PHGEOMETRY_PHGEOMUTILITY_H
#define PHGEOMETRY_PHGEOMUTILITY_H

#include <string>

class PHCompositeNode;
class TGeoManager;
class PHGeomTGeo;
class PHGeomIOTGeo;

//! Toolsets to do geometry operations
class PHGeomUtility
{
 public:
  //! Main user interface: DST node -> TGeoManager for downstream use
  static TGeoManager *
  GetTGeoManager(PHCompositeNode *topNode);

  //! TGeo ROOT/GDML/Macro file -> DST node with automatic file type discrimination based on file names
  static int
  ImportGeomFile(PHCompositeNode *topNode, const std::string &geometry_file);

  //! gGeoManager -> DST node
  static int
  ImportCurrentTGeoManager(PHCompositeNode *topNode);

  //! DST node -> TGeoManager -> export files, like gdml, .root or .C formats
  static void
  ExportGeomtry(PHCompositeNode *topNode, const std::string &geometry_file);

  //! Get non-persistent PHGeomTGeo from DST nodes. If not found, make a new one
  static PHGeomTGeo *
  GetGeomTGeoNode(PHCompositeNode *topNode, bool build_new = true);

  //! Get persistent PHGeomIOTGeo from DST nodes. If not found, make a new one
  static PHGeomIOTGeo *
  GetGeomIOTGeoNode(PHCompositeNode *topNode, bool build_new = true);

  //! Update persistent PHGeomIOTGeo node RUN/GEOMETRY_IO based on run-time object PHGeomTGeo at RUN/GEOMETRY
  //! \return the updated PHGeomIOTGeo from DST tree
  static PHGeomIOTGeo *
  UpdateIONode(PHCompositeNode *topNode);

  //! Build or update PHGeomTGeo node RUN/GEOMETRY based on the persistent PHGeomIOTGeo node RUN/GEOMETRY_IO
  //! \return the updated PHGeomTGeo from DST tree
  static PHGeomTGeo *
  LoadFromIONode(PHCompositeNode *topNode);

  //! Make a name for tmp geometry file
  //! Geometry files gain a size of ~10MB and it used in translation from Geant4 to DST format.
  //! This tmp file should be on a local file system (/tmp/) and write/deletable
  static std::string
  GenerateGeometryFileName(const std::string &filename_extension = "gdml");

  //! delete the geometry file after use
  static bool
  RemoveGeometryFile(const std::string &file_name);

  //! Verbosity for geometry IO like, TGeoMangers
  static void SetVerbosity(int v);

  //! Verbosity for geometry IO like, TGeoMangers
  static int GetVerbosity();

  //! DST node name for RunTime geometry object
  static std::string
  GetDSTNodeName()
  {
    return std::string("GEOMETRY");
  }

  //! DST node name for persistent geometry storage node
  static std::string
  GetDSTIONodeName()
  {
    return std::string("GEOMETRY_IO");
  }

  //! Base path name for temp geometry GDML file used in GenerateGeometryFileName().
  //! User can overwrite it to e.g. local directory with  PHGeomUtility::SetGenerateGeometryFileNameBase('./');
  static void SetGenerateGeometryFileNameBase(const std::string &base) { mg_GenerateGeometryFileNameBase = base; }

 private:
  PHGeomUtility() = delete;
  ~PHGeomUtility() = delete;

  //! Base path name for temp geometry GDML file used in GenerateGeometryFileName().
  //! User can overwrite it to e.g. local directory with  PHGeomUtility::SetGenerateGeometryFileNameBase('./');
  static std::string mg_GenerateGeometryFileNameBase;
};

#endif
