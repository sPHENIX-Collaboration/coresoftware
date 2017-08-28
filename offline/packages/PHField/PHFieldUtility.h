#ifndef PHFieldUtility_HH__
#define PHFieldUtility_HH__

#include <ctime>
#include <iostream>
#include <map>
#include <set>
#include <string>

#include "PHFieldConfig.h"

class PHCompositeNode;
class PHField;
class PHFieldConfig;

//! Toolsets to do geometry operations
class PHFieldUtility
{
 protected:
  // static tool sets only
  PHFieldUtility();
  virtual ~PHFieldUtility();

 public:

  //! Make a default PHFieldConfig
  static PHFieldConfig * MakeDefaultFieldConfig();

  //! Get transient PHField from DST nodes. If not found, make a new one based on default_config
  static PHField *
  GetFieldMapNode(PHFieldConfig *default_config = nullptr, PHCompositeNode *topNode = nullptr);

  //! Get persistent PHFieldConfig from DST nodes. If not found, make a new one based on default_config
  static PHFieldConfig *
  GetGeomIOTGeoNode(PHFieldConfig *default_config = nullptr, PHCompositeNode *topNode = nullptr);

  //! Build or update PHField node RUN/GEOMETRY based on the persistent PHGeomIOTGeo node RUN/GEOMETRY_IO
  //! \return the updated PHGeomTGeo from DST tree
  static PHField *
  LoadFromIONode(PHCompositeNode *topNode);

  //! DST node name for RunTime field map object
  static std::string
  GetDSTFieldMapNodeName()
  {
    return std::string("FIELD_MAP");
  }

  //! DST node name for persistent field configuration node
  static std::string
  GetDSTConfigNodeName()
  {
    return std::string("FIELD_CONFIG");
  }

 protected:

};

#endif /* PHFieldUtility_HH__ */
