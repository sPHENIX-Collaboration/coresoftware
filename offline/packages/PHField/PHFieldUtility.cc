#include "PHFieldUtility.h"

#include "PHField.h"
#include "PHField2D.h"
#include "PHField3DCartesian.h"
#include "PHField3DCylindrical.h"
#include "PHFieldBeast.h"
#include "PHFieldCleo.h"
#include "PHFieldConfig.h"
#include "PHFieldConfigv1.h"
#include "PHFieldUniform.h"

#include <fun4all/Fun4AllServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TSystem.h>

#include <cassert>
#include <cstdlib>  // for getenv
#include <iostream>

using namespace std;

PHField *
PHFieldUtility::BuildFieldMap(const PHFieldConfig *field_config, const int verbosity)
{
  assert(field_config);

  if (verbosity)
  {
    cout << "PHFieldUtility::BuildFieldMap - construction field with configuration: ";
    field_config->identify();
  }

  PHField *field(nullptr);

  switch (field_config->get_field_config())
  {
  case PHFieldConfig::kFieldUniform:
    //    return "Constant field";

    field = new PHFieldUniform(
        field_config->get_field_mag_x(),
        field_config->get_field_mag_y(),
        field_config->get_field_mag_z());

    break;
  case PHFieldConfig::kField2D:
    //    return "2D field map expressed in cylindrical coordinates";
    field = new PHField2D(
        field_config->get_filename(),
        verbosity,
        field_config->get_magfield_rescale());
    break;

  case PHFieldConfig::kField3DCylindrical:
    //    return "3D field map expressed in cylindrical coordinates";
    field = new PHField3DCylindrical(
        field_config->get_filename(),
        verbosity,
        field_config->get_magfield_rescale());
    break;

  case PHFieldConfig::Field3DCartesian:
    //    return "3D field map expressed in Cartesian coordinates";
    field = new PHField3DCartesian(
        field_config->get_filename(),
        field_config->get_magfield_rescale());
    break;

  case PHFieldConfig::kFieldBeast:
    //    return "2D Beast field map expressed in Cartesian coordinates";
    cout << "calling PHFieldBeast scale " << field_config->get_magfield_rescale() << endl;
    field = new PHFieldBeast(
        field_config->get_filename(),
        verbosity,
        field_config->get_magfield_rescale());
    break;

  case PHFieldConfig::kFieldCleo:
    //    return "2D Beast field map expressed in Cartesian coordinates";
    cout << "calling PHFieldCleo scale " << field_config->get_magfield_rescale() << endl;
    field = new PHFieldCleo(
        field_config->get_filename(),
        verbosity,
        field_config->get_magfield_rescale());
    break;

  default:
    cout << "PHFieldUtility::BuildFieldMap - Invalid Field Configuration" << endl;
    //    return nullptr;
    //    return "Invalid Field";
  }
  assert(field);  // Check for Invalid Field
  return field;
}

//! Make a default PHFieldConfig
//! Field map = /phenix/upgrades/decadal/fieldmaps/sPHENIX.2d.root
//! Field Scale to 1.4/1.5
//! \output default field configuration object. Caller assumes ownership
PHFieldConfig *PHFieldUtility::DefaultFieldConfig()
{
  return new PHFieldConfigv1(PHFieldConfigv1::kField2D,
                             (string(getenv("CALIBRATIONROOT")) + string("/Field/Map/sPHENIX.2d.root")),
                             -1.4 / 1.5);
}

//! Get transient PHField from DST nodes. If not found, make a new one based on default_config
PHField *
PHFieldUtility::GetFieldMapNode(const PHFieldConfig *default_config, PHCompositeNode *topNode, const int verbosity)
{
  if (topNode == nullptr)
  {
    topNode = Fun4AllServer::instance()->topNode();
  }
  assert(topNode);
  PHNodeIterator iter(topNode);

  // Looking for the RUN node
  PHCompositeNode *parNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  if (!parNode)
  {
    cout << PHWHERE << ": PAR Node missing, request aborting.";
    gSystem->Exit(1);
  }

  PHField *field = findNode::getClass<PHField>(parNode, GetDSTFieldMapNodeName());
  if (!field)
  {
    PHFieldConfig *field_config = GetFieldConfigNode(default_config, topNode, verbosity);
    assert(field_config);

    field = BuildFieldMap(field_config, verbosity > 0 ? verbosity - 1 : verbosity);
    assert(field);

    parNode->addNode(new PHDataNode<PHField>(field, GetDSTFieldMapNodeName()));
  }

  return field;
}

//! Get persistent PHGeomIOTGeo from DST nodes. If not found, make a new one
PHFieldConfig *
PHFieldUtility::GetFieldConfigNode(const PHFieldConfig *default_config, PHCompositeNode *topNode, const int verbosity)
{
  if (topNode == nullptr)
  {
    topNode = Fun4AllServer::instance()->topNode();
  }
  assert(topNode);

  PHNodeIterator iter(topNode);

  // Looking for the RUN node
  PHCompositeNode *runNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    cout << PHWHERE << ": RUN Node missing, aborting.";
    gSystem->Exit(1);
  }

  PHFieldConfig *field = findNode::getClass<PHFieldConfig>(runNode, GetDSTConfigNodeName());
  if (!field)
  {
    if (!default_config)
    {
      field = DefaultFieldConfig();
      if (verbosity)
      {
        cout << "PHFieldUtility::GetFieldConfigNode - field map with configuration from build-in default: ";
        field->identify();
      }
    }
    else
    {
      field = static_cast<PHFieldConfig *>(default_config->CloneMe());
      if (verbosity)
      {
        cout << "PHFieldUtility::GetFieldConfigNode - field map with configuration from input default: ";
        field->identify();
      }
    }

    assert(field);
    runNode->addNode(new PHIODataNode<PHObject>(field, GetDSTConfigNodeName(), "PHObject"));
  }
  else
  {
    if (verbosity)
    {
      cout << "PHFieldUtility::GetFieldConfigNode - field map with configuration from DST/RUN: ";
      field->identify();
    }
  }

  return field;
}
