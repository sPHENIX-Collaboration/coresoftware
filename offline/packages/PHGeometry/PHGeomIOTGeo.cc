// $Id: $                                                                                             

/*!
 * \file PHGeomIOTGeo.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHGeomIOTGeo.h"

#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TMemFile.h>

#include <cassert>
#include <sstream>
#include <iostream>

using namespace std;

PHGeomIOTGeo::PHGeomIOTGeo() :
    Data(0)
{
}

PHGeomIOTGeo::~PHGeomIOTGeo()
{
  Reset();
}

PHObject*
PHGeomIOTGeo::clone() const
{
  PHGeomIOTGeo * geo = new PHGeomIOTGeo(*this);
  return geo;
}

void
PHGeomIOTGeo::SetGeometry(const TGeoVolume * g)
{
  if (!g)
    {
      cout << __PRETTY_FUNCTION__ << " - Error - Invalid input" << endl;
      return;
    }

  // Stream TGeoVolume into binary stream with its streamer using TFIle utility
  TMemFile f1("mem","CREATE");
  g->Write("TOP");
  f1.Close();

  const Long64_t n = f1.GetSize();

  Data.resize(n);
  Long64_t n1 = f1.CopyTo(Data.data(), n);
  assert(n1 == n);

}

TGeoVolume *
PHGeomIOTGeo::GetGeometryCopy()
{
  if (not isValid()) return nullptr;

  TMemFile f2("mem2", Data.data(), Data.size(), "READ");
  TGeoVolume * vol = dynamic_cast<TGeoVolume *>(f2.Get("TOP"));
  assert(vol);
  f2.Close();

  return vol;
}

TGeoManager *
PHGeomIOTGeo::
ConstructTGeoManager()
{
  if (not isValid()) return nullptr;

  // build new TGeoManager
  TGeoManager * tgeo = new TGeoManager("PHGeometry", "");
  assert(tgeo);

  TGeoVolume * vol = GetGeometryCopy();
  vol->RegisterYourself();

  tgeo->SetTopVolume(vol);
//  tgeo->CloseGeometry();

  ostringstream stitle;
  stitle
      << "TGeoManager built by PHGeomUtility::LoadFromIONode based on RUN/GEOMETRY_IO node with name ("
      << vol->GetName() << ") and title ("
      << vol->GetTitle() << ")";

  tgeo->SetTitle(stitle.str().c_str());

  return tgeo;
}

/** identify Function from PHObject
 @param os Output Stream
 */
void
PHGeomIOTGeo::identify(std::ostream& os) const
{
  os << "PHGeomIOTGeo - ";
  if (isValid())
    os << " with geometry data " << Data.size()<<"Byte";
  else
    os << "Empty";
  os << endl;
}

/// Clear Event
void
PHGeomIOTGeo::Reset()
{
  Data.resize(0);
}

/// isValid returns non zero if object contains vailid data
int
PHGeomIOTGeo::isValid() const
{
  return Data.size();
}
