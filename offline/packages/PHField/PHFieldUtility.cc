#include "PHFieldUtility.h"

#include "PHField.h"
#include "PHField2D.h"
#include "PHField3DCartesian.h"
#include "PHField3DCylindrical.h"
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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <TSystem.h>
#pragma GCC diagnostic pop

#include <cassert>
#include <cstdlib>  // for getenv
#include <iostream>

PHField *
PHFieldUtility::BuildFieldMap(const PHFieldConfig *field_config, float inner_radius, float outer_radius, float size_z, const int verbosity)
{
  assert(field_config);

  if (verbosity)
  {
    std::cout << "PHFieldUtility::BuildFieldMap - construction field with configuration: ";
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
        field_config->get_magfield_rescale(),
	inner_radius,
	outer_radius,
        size_z);
    break;

  default:
    std::cout << "PHFieldUtility::BuildFieldMap - Invalid Field Configuration: " << field_config->get_field_config() << std::endl;
    gSystem->Exit(1);
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
  char *calibrationroot = getenv("CALIBRATIONROOT");
  std::string fieldmap = "sphenix3dbigmapxyz_gap_rebuild.root";
  if (calibrationroot != nullptr)
  {
    fieldmap = std::string(calibrationroot) + "/Field/Map/" + fieldmap;
  }
  return new PHFieldConfigv1(PHFieldConfigv1::Field3DCartesian, fieldmap, 1.);
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
    std::cout << PHWHERE << ": PAR Node missing, request aborting.";
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
    std::cout << PHWHERE << ": RUN Node missing, aborting.";
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
        std::cout << "PHFieldUtility::GetFieldConfigNode - field map with configuration from build-in default: ";
        field->identify();
      }
    }
    else
    {
      field = static_cast<PHFieldConfig *>(default_config->CloneMe());
      if (verbosity)
      {
        std::cout << "PHFieldUtility::GetFieldConfigNode - field map with configuration from input default: ";
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
      std::cout << "PHFieldUtility::GetFieldConfigNode - field map with configuration from DST/RUN: ";
      field->identify();
    }
  }

  return field;
}
