#include "PHGeomUtility.h"
#include "PHGeomTGeo.h"

#include <TGeoManager.h>
#include <TROOT.h>

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>

#include <sys/types.h>
#include <unistd.h> // for generate unique local file
using namespace std;

PHGeomUtility::PHGeomUtility()
{
}

PHGeomUtility::~PHGeomUtility()
{

}

//! DST node -> TGeoManager for downstream use
TGeoManager *
PHGeomUtility::GetTGeoManager(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the RUN node
  PHCompositeNode *runNode = static_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  if (!runNode)
    {
      cout << __PRETTY_FUNCTION__ << " - Error - RUN Node missing." << endl;
      return NULL;
    }

  PHGeomTGeo *dst_geom = findNode::getClass<PHGeomTGeo>(runNode,
      GetDSTNodeName());
  if (!dst_geom)
    {
      cout << __PRETTY_FUNCTION__ << " - Error - RUN Geometry Node missing."
          << endl;
      return NULL;
    }

  return dst_geom->GetGeometry();
}

int
PHGeomUtility::ImportGeomFile(PHCompositeNode *topNode,
    const std::string & geometry_file)
{

  PHNodeIterator iter(topNode);

  // Looking for the RUN node
  PHCompositeNode *runNode = static_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  if (!runNode)
    {
      cout << __PRETTY_FUNCTION__ << "RUN Node missing, request aborting."
          << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  PHGeomTGeo *dst_geom = findNode::getClass<PHGeomTGeo>(runNode,
      GetDSTNodeName());
  if (!dst_geom)
    {
      dst_geom = new PHGeomTGeo();
      PHDataNode<PHObject> *GeomNode = new PHDataNode<PHObject>(dst_geom,
          GetDSTNodeName(), "PHObject");
      runNode->addNode(GeomNode);
    }

  dst_geom->SetGeometry(TGeoManager::Import(geometry_file.c_str()));

  if (dst_geom->GetGeometry() == NULL)
    {
      cout << __PRETTY_FUNCTION__ << "failed to import " << geometry_file
          << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

//int
//PHGeomUtility::ImportGDML(PHCompositeNode *topNode, const std::string & gdml_file)
//{
//  PHNodeIterator iter(topNode);
//
//  // Looking for the RUN node
//  PHCompositeNode *runNode = static_cast<PHCompositeNode*>(iter.findFirst(
//      "PHCompositeNode", "RUN"));
//  if (!runNode)
//    {
//      cout << __PRETTY_FUNCTION__ << "RUN Node missing, request aborting."
//          << endl;
//      return Fun4AllReturnCodes::ABORTRUN;
//    }
//
////  runNode->
//  PHGeomTGeo *dst_geom = findNode::getClass<PHGeomTGeo>(runNode,
//      GetDSTNodeName());
//  if (!dst_geom)
//    {
//      dst_geom = new PHGeomTGeo();
//      PHIODataNode<PHObject> *GeomNode = new PHIODataNode<PHObject>(dst_geom,
//          GetDSTNodeName(), "PHObject");
//      runNode->addNode(GeomNode);
//    }
//
//  //TODO: GDML import
////  dst_geom->SetGeometry(TGeoManager::Import(gdml_file.c_str()));
//  if (gGeoManager) delete gGeoManager;
//
//
//
//  if (dst_geom->GetGeometry() == NULL)
//    {
//      cout << __PRETTY_FUNCTION__ << "failed to import " << gdml_file
//          << endl;
//      return Fun4AllReturnCodes::ABORTRUN;
//    }
//
//
//  return Fun4AllReturnCodes::EVENT_OK;
//}

std::string
PHGeomUtility::GenerateGeometryFileName(const std::string & filename_extension)
{
  stringstream file;
  file << "/tmp/" << "PHGeomUtility_geom_file_" << ::getpid() << "."
      << filename_extension << endl;

  return file.str();
}

//! delete the geometry file after use
bool
PHGeomUtility::RemoveGeometryFile(const std::string & file_name)
{
  ifstream ifile(file_name);

  if (ifile)
    {
      ifile.close();
      if (remove(file_name.c_str()) != 0)
        {
          cout << __PRETTY_FUNCTION__ << " - Error - can not remove file "
              << file_name << endl;
          return false;
        }
      else
        return true;
    }
  else
    return true; // file do not exist

}
