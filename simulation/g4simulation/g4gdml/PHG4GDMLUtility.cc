// $Id: $

/*!
 * \file PHG4GDMLUtility.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4GDMLUtility.hh"
#include "PHG4GDMLWriteStructure.hh"

#include <Geant4/G4VPhysicalVolume.hh>

#include <cassert>

PHG4GDMLUtility::~PHG4GDMLUtility()
{
  // TODO Auto-generated destructor stub
}
PHG4GDMLUtility::PHG4GDMLUtility()
{
  // TODO Auto-generated constructor stub
}

void PHG4GDMLUtility::Dump_GDML(const std::string &filename, G4VPhysicalVolume * vol)
{
  PHG4GDMLWriteStructure gdml_parser;
  assert(vol);
  assert(vol->GetLogicalVolume());
  gdml_parser.Write(filename, vol->GetLogicalVolume(),get_PHG4GDML_Schema(), 0, true);
}
