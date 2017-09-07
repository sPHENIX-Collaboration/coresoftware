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
  //! Make a default PHFieldConfig as in default macro of pro.3 release
  //! Field map = /phenix/upgrades/decadal/fieldmaps/sPHENIX.2d.root
  //! Field Scale to 1.4/1.5
  //! \output default field configuration object. Caller assumes ownership
  static PHFieldConfig *DefaultFieldConfig();

  //! Get transient PHField from DST nodes. If not found, make a new one based on default_config
  //! \param[in]  default_config  default configuraiton if not on DST. If nullptr, use DefaultFieldConfig() as the default
  //! \param[in]  topNode         you know who....
  static PHField *
  GetFieldMapNode(const PHFieldConfig *default_config = nullptr, PHCompositeNode *topNode = nullptr, const int verbosity = 0);

  //! Get persistent PHFieldConfig from DST nodes. If not found, make a new one based on default_config
  //! \param[in]  default_config  default configuraiton if not on DST. If nullptr, use DefaultFieldConfig() as the default
  //! \param[in]  topNode         you know who....
  static PHFieldConfig *
  GetFieldConfigNode(const PHFieldConfig *default_config = nullptr, PHCompositeNode *topNode = nullptr, const int verbosity = 0);

  //! Build or build field map with a configuration object
  static PHField *
  BuildFieldMap(const PHFieldConfig *field_config, const int verbosity = 0);

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
