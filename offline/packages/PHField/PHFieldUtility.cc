#include "PHFieldUtility.h"

#include "PHFieldConfig_v1.h"

#include "PHField2D.h"
#include "PHField3DCartesian.h"
#include "PHField3DCylindrical.h"
#include "PHFieldConst.h"

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <sys/types.h>
#include <unistd.h>  // for generate unique local file
using namespace std;

PHFieldUtility::PHFieldUtility()
{
}

PHFieldUtility::~PHFieldUtility()
{
}

PHField *
PHFieldUtility::BuildFieldMap(const PHFieldConfig *field_config)
{
  assert(field_config);

  PHField *field(nullptr);

  switch (get_field_config())
  {
  case kFieldConstant:
    //    return "Constant field";
    break;
  case kField2D:
    //    return "2D field map expressed in cylindrical coordinates";
    break;
  case kField3DCylindrical:
    //    return "3D field map expressed in cylindrical coordinates";
    break;
  case Field3DCartesian:
    //    return "3D field map expressed in Cartesian coordinates";

    break;
  default:
    cout << "PHFieldUtility::BuildFieldMap - Invalid Field Configuration" << endl;
    assert(0);  // Invalid Field
                //    return nullptr;
                //    return "Invalid Field";
  }
}

//! Make a default PHFieldConfig
//! Field map = /phenix/upgrades/decadal/fieldmaps/sPHENIX.2d.root
//! Field Scale to 1.4/1.5
//! \output default field configuration object. Caller assumes ownership
PHFieldConfig *
PHFieldUtility::
    DefaultFieldConfig()
{
  return new PHFieldConfig_v1(PHFieldConfig_v1::kField2D,
                              "/phenix/upgrades/decadal/fieldmaps/sPHENIX.2d.root",
                              1.4 / 1.5);
}

//! Get transient PHField from DST nodes. If not found, make a new one based on default_config
static PHField *
PHFieldUtility::GetFieldMapNode(PHFieldConfig *default_config, PHCompositeNode *topNode)
{
  if (topNode == nullptr) topNode = Fun4AllServer::instance()->topNode();
  assert(topNode);
  PHNodeIterator iter(topNode);

  // Looking for the RUN node
  PHCompositeNode *parNode = static_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "PAR"));
  if (!parNode)
  {
    stringstream serr;
    serr << __PRETTY_FUNCTION__ << ": PAR Node missing, request aborting.";
    cout << serr.str() << endl;

    throw runtime_error(serr.str());

    return NULL;
  }

  PHField *field = findNode::getClass<PHField>(parNode, GetDSTFieldMapNodeName());
  if (!field and default_config)
  {
    PHFieldConfig *field_config =
        GetFieldConfigNode(default_config, topNode);
    assert(field_config);

    field = BuildFieldMap(field_config);
    assert(field);

    parNode->addNode(new PHDataNode<PHObject>(field, GetDSTFieldMapNodeName(), "PHObject"));
  }

  return field;
}

//! Get persistent PHGeomIOTGeo from DST nodes. If not found, make a new one
PHFieldConfig *
PHFieldUtility::GetFieldConfigNode(PHFieldConfig *default_config, PHCompositeNode *topNode)
{
  if (topNode == nullptr) topNode = Fun4AllServer::instance()->topNode();
  assert(topNode);

  PHNodeIterator iter(topNode);

  // Looking for the RUN node
  PHCompositeNode *runNode = static_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  if (!runNode)
  {
    stringstream serr;
    serr << __PRETTY_FUNCTION__ << ": RUN Node missing, request aborting.";
    cout << serr.str() << endl;

    throw runtime_error(serr.str());

    return NULL;
  }

  PHField *field = findNode::getClass<PHField>(parNode,
                                               GetDSTConfigNodeName());
  if (!field )
  {
    if (!default_config)
      field = DefaultFieldConfig();
    else
      field = static_cast<PHFieldConfig *>(default_config->clone());

    assert(field);
    runNode->addNode(new PHDataNode<PHObject>(field,
                                              GetDSTConfigNodeName(), "PHObject"));
  }

  return field;
}




