#include "PHGeomUtility.h"
#include "PHGeomTGeo.h"

#include <TGeoManager.h>

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <cassert>
#include <iostream>

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
      cout << __PRETTY_FUNCTION__ << " - Error - RUN Node missing."
          << endl;
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
    const std::string & geometry_root_file)
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
      PHIODataNode<PHObject> *GeomNode = new PHIODataNode<PHObject>(dst_geom,
          GetDSTNodeName(), "PHObject");
      runNode->addNode(GeomNode);
    }

  dst_geom->SetGeometry(TGeoManager::Import(geometry_root_file.c_str()));

  if (dst_geom->GetGeometry() == NULL)
    {
      cout << __PRETTY_FUNCTION__ << "failed to import " << geometry_root_file
          << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int
PHGeomUtility::ImportGDML(PHCompositeNode *topNode, const std::string & gdml_file)
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

//  runNode->
  PHGeomTGeo *dst_geom = findNode::getClass<PHGeomTGeo>(runNode,
      GetDSTNodeName());
  if (!dst_geom)
    {
      dst_geom = new PHGeomTGeo();
      PHIODataNode<PHObject> *GeomNode = new PHIODataNode<PHObject>(dst_geom,
          GetDSTNodeName(), "PHObject");
      runNode->addNode(GeomNode);
    }

  //TODO: GDML import
  dst_geom->SetGeometry(TGeoManager::Import(gdml_file.c_str()));

  if (dst_geom->GetGeometry() == NULL)
    {
      cout << __PRETTY_FUNCTION__ << "failed to import " << gdml_file
          << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }


  return Fun4AllReturnCodes::EVENT_OK;
}
