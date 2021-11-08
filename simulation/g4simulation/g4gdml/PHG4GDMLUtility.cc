// $Id: $

/*!
 * \file PHG4GDMLUtility.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4GDMLUtility.hh"
#include "PHG4GDMLConfig.hh"
#include "PHG4GDMLWriteStructure.hh"

#include <fun4all/Fun4AllServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>

#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4GDMLWriteStructure.hh>

#include <cassert>
#include <iostream>  // for operator<<, stringstream
#include <sstream>
#include <stdexcept>

using namespace std;

void PHG4GDMLUtility::Dump_GDML(const std::string &filename, G4VPhysicalVolume *vol, PHCompositeNode *topNode)
{
  if (topNode == nullptr)
  {
    Fun4AllServer *se = Fun4AllServer::instance();
    topNode = se->topNode();
  }

  const PHG4GDMLConfig *config =
      GetOrMakeConfigNode(topNode);
  assert(config);

  PHG4GDMLWriteStructure gdml_parser(config);
  assert(vol);
  assert(vol->GetLogicalVolume());

  xercesc::XMLPlatformUtils::Initialize();
  gdml_parser.Write(filename, vol->GetLogicalVolume(), get_PHG4GDML_Schema(), 0, true);
  xercesc::XMLPlatformUtils::Terminate();
}

void PHG4GDMLUtility::Dump_G4_GDML(const std::string &filename, G4VPhysicalVolume *vol)
{

  G4GDMLWriteStructure gdml_parser;
  assert(vol);
  assert(vol->GetLogicalVolume());

  xercesc::XMLPlatformUtils::Initialize();
  gdml_parser.Write(filename, vol->GetLogicalVolume(), get_PHG4GDML_Schema(), 0, true);
  xercesc::XMLPlatformUtils::Terminate();
}

PHG4GDMLConfig *PHG4GDMLUtility::GetOrMakeConfigNode(PHCompositeNode *topNode, bool build_new)
{
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

    return nullptr;
  }

  PHG4GDMLConfig *config = findNode::getClass<PHG4GDMLConfig>(parNode,
                                                              getDSTNodeName());
  if (!config and build_new)
  {
    config = new PHG4GDMLConfig();
    PHDataNode<PHObject> *node = new PHDataNode<PHObject>(config,
                                                          getDSTNodeName(), "PHObject");
    parNode->addNode(node);
  }

  return config;
}
