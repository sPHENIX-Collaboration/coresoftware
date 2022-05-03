#include "PHGeomUtility.h"

#include "PHGeomIOTGeo.h"
#include "PHGeomTGeo.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <TGeoManager.h>

#include <uuid/uuid.h>

#include <unistd.h>  // for generate unique local file
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//! DST node -> TGeoManager for downstream use
TGeoManager *
PHGeomUtility::GetTGeoManager(PHCompositeNode *topNode)
{
  PHGeomTGeo *dst_geom = GetGeomTGeoNode(topNode, true);
  if (!dst_geom)
  {
    cout << __PRETTY_FUNCTION__
         << " - Error - Can NOT construct geometry node." << endl;
    exit(1);
    return nullptr;
  }

  if (not dst_geom->isValid())
  {
    // try to construct the geometry node
    dst_geom = LoadFromIONode(topNode);
  }

  if (TGeoManager::GetDefaultUnits() != TGeoManager::kRootUnits )
  {
    cout << __PRETTY_FUNCTION__ << " TGeoManager was not constructed with RootUnits, which potentially leads to unit mismatch with Fun4All. This is considered a fatal error."
        <<endl;

    exit(1);
    return nullptr;
  }

  UpdateIONode(topNode);

  return dst_geom->GetGeometry();
}

void PHGeomUtility::ExportGeomtry(PHCompositeNode *topNode,
                                  const std::string &geometry_file)
{
  TGeoManager *tgeo = GetTGeoManager(topNode);

  assert(tgeo);

  tgeo->Export(geometry_file.c_str());
}

int PHGeomUtility::ImportGeomFile(PHCompositeNode *topNode,
                                  const std::string &geometry_file)
{
  PHGeomTGeo *dst_geom = GetGeomTGeoNode(topNode);
  assert(dst_geom);

  dst_geom->Reset();

  TGeoManager::SetVerboseLevel(GetVerbosity());

  // force TGeoManager to use the Fun4All unit of cm
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,23,2)
  TGeoManager::LockDefaultUnits(kFALSE);
  TGeoManager::SetDefaultUnits( TGeoManager::kRootUnits );
  TGeoManager::LockDefaultUnits(kTRUE);
#else
  TGeoManager::SetDefaultRootUnits();
#endif

  dst_geom->SetGeometry(TGeoManager::Import(geometry_file.c_str()));

  if (dst_geom->GetGeometry() == nullptr)
  {
    cout << __PRETTY_FUNCTION__ << "failed to import " << geometry_file
         << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  UpdateIONode(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHGeomUtility::ImportCurrentTGeoManager(PHCompositeNode *topNode)
{
  PHGeomTGeo *dst_geom = GetGeomTGeoNode(topNode);
  assert(dst_geom);

  if (dst_geom->GetGeometry() == gGeoManager)
    return Fun4AllReturnCodes::EVENT_OK;  // noting to be done

  assert(not dst_geom->isValid());  // check that it is uninitialized
  dst_geom->SetGeometry(gGeoManager);
  TGeoManager::SetVerboseLevel(GetVerbosity());

  return Fun4AllReturnCodes::EVENT_OK;
}

//! Get non-persistent PHGeomTGeo from DST nodes. If not found, make a new one
PHGeomTGeo *
PHGeomUtility::GetGeomTGeoNode(PHCompositeNode *topNode, bool build_new)
{
  PHNodeIterator iter(topNode);

  // Looking for the RUN node
  PHCompositeNode *parNode = static_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "PAR"));
  if (!parNode)
  {
    ostringstream serr;
    serr << __PRETTY_FUNCTION__ << ": PAR Node missing, request aborting.";
    cout << serr.str() << endl;

    throw runtime_error(serr.str());

    return nullptr;
  }

  PHGeomTGeo *dst_geom = findNode::getClass<PHGeomTGeo>(parNode,
                                                        GetDSTNodeName());
  if (!dst_geom and build_new)
  {
    dst_geom = new PHGeomTGeo();
    PHDataNode<PHObject> *GeomNode = new PHDataNode<PHObject>(dst_geom,
                                                              GetDSTNodeName(), "PHObject");
    parNode->addNode(GeomNode);
  }

  return dst_geom;
}

//! Get persistent PHGeomIOTGeo from DST nodes. If not found, make a new one
PHGeomIOTGeo *
PHGeomUtility::GetGeomIOTGeoNode(PHCompositeNode *topNode, bool build_new)
{
  PHNodeIterator iter(topNode);

  // Looking for the RUN node
  PHCompositeNode *runNode = static_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  if (!runNode)
  {
    ostringstream serr;
    serr << __PRETTY_FUNCTION__ << ": RUN Node missing, request aborting.";
    cout << serr.str() << endl;

    throw runtime_error(serr.str());

    return nullptr;
  }

  PHGeomIOTGeo *dst_geom = findNode::getClass<PHGeomIOTGeo>(runNode,
                                                            GetDSTIONodeName());
  if (!dst_geom and build_new)
  {
    dst_geom = new PHGeomIOTGeo();
    PHIODataNode<PHObject> *GeomNode = new PHIODataNode<PHObject>(dst_geom,
                                                                  GetDSTIONodeName(), "PHObject");
    runNode->addNode(GeomNode);
  }

  return dst_geom;
}

std::string
PHGeomUtility::GenerateGeometryFileName(const std::string &filename_extension)
{
  ostringstream file;

  uuid_t uu;
  uuid_generate(uu);
  char uuid[50];
  uuid_unparse(uu, uuid);

  file << mg_GenerateGeometryFileNameBase << "/"
       << "PHGeomUtility_geom_file_" << uuid << "."
       << filename_extension;

  return file.str();
}

std::string PHGeomUtility::mg_GenerateGeometryFileNameBase = "/tmp";

//! delete the geometry file after use
bool PHGeomUtility::RemoveGeometryFile(const std::string &file_name)
{
  fstream ifile(file_name, ios_base::in);

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
    return true;  // file do not exist
}

//! Update persistent PHGeomIOTGeo node RUN/GEOMETRY_IO based on run-time object PHGeomTGeo at RUN/GEOMETRY
//! \return the updated PHGeomIOTGeo from DST tree
PHGeomIOTGeo *
PHGeomUtility::UpdateIONode(PHCompositeNode *topNode)
{
  PHGeomTGeo *dst_geom = GetGeomTGeoNode(topNode, false);

  if (not dst_geom)
  {
    cout << __PRETTY_FUNCTION__
         << " - ERROR - failed to update PHGeomIOTGeo node RUN/GEOMETRY_IO due to missing PHGeomTGeo node at RUN/GEOMETRY"
         << endl;
    return nullptr;
  }
  if (not dst_geom->isValid())
  {
    cout << __PRETTY_FUNCTION__
         << " - ERROR - failed to update PHGeomIOTGeo node RUN/GEOMETRY_IO due to invalid PHGeomTGeo node at RUN/GEOMETRY"
         << endl;
    return nullptr;
  }

  PHGeomIOTGeo *dst_geom_io = GetGeomIOTGeoNode(topNode, true);
  assert(dst_geom_io);

  dst_geom_io->SetGeometry(dst_geom->GetGeometry()->GetTopVolume());

  return dst_geom_io;
}

//! Build or update PHGeomTGeo node RUN/GEOMETRY based on the persistent PHGeomIOTGeo node RUN/GEOMETRY_IO
//! \return the updated PHGeomTGeo from DST tree
PHGeomTGeo *
PHGeomUtility::LoadFromIONode(PHCompositeNode *topNode)
{
  PHGeomIOTGeo *dst_geom_io = GetGeomIOTGeoNode(topNode, false);

  if (not dst_geom_io)
  {
    cout << __PRETTY_FUNCTION__
         << " - ERROR - failed to update PHGeomTGeo node RUN/GEOMETRY due to missing PHGeomIOTGeo node at RUN/GEOMETRY_IO"
         << endl;
    return nullptr;
  }
  if (not dst_geom_io->isValid())
  {
    cout << __PRETTY_FUNCTION__
         << " - ERROR - failed to update PHGeomTGeo node RUN/GEOMETRY due to invalid PHGeomIOTGeo node at RUN/GEOMETRY_IO"
         << endl;
    return nullptr;
  }

  // build new TGeoManager
  TGeoManager::SetVerboseLevel(GetVerbosity());
  TGeoManager *tgeo = dst_geom_io->ConstructTGeoManager();
  tgeo->CloseGeometry();

  PHGeomTGeo *dst_geom = GetGeomTGeoNode(topNode, true);
  assert(dst_geom);
  dst_geom->SetGeometry(tgeo);

  return dst_geom;
}

//! Verbosity for geometry IO like, TGeoMangers
void PHGeomUtility::SetVerbosity(int v)
{
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("PHGEOMETRY_VERBOSITY", v);
}

//! Verbosity for geometry IO like, TGeoMangers
int PHGeomUtility::GetVerbosity()
{
  recoConsts *rc = recoConsts::instance();
  if (rc->FlagExist("PHGEOMETRY_VERBOSITY"))
    return rc->get_IntFlag("PHGEOMETRY_VERBOSITY");
  else
    return 0;
}
