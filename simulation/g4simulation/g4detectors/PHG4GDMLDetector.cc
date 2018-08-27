// $Id: $

/*!
 * \file PHG4GDMLDetector.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4GDMLDetector.h"

#include <phparameter/PHParameters.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4SystemOfUnits.hh>

using namespace std;

PHG4GDMLDetector::PHG4GDMLDetector(PHCompositeNode *Node, const std::string &dnam, PHParameters *parameters)
  : PHG4Detector(Node, dnam)
  , m_GDMPath(parameters->get_string_param("GDMPath"))
  , m_TopVolName(parameters->get_string_param("TopVolName"))
  , m_placeX(parameters->get_double_param("place_y") * cm)
  , m_placeY(parameters->get_double_param("inner_radius") * cm)
  , m_placeZ(parameters->get_double_param("place_z") * cm)
  , m_rotationX(parameters->get_double_param("rot_x") * rad)
  , m_rotationY(parameters->get_double_param("rot_y") * rad)
  , m_rotationZ(parameters->get_double_param("rot_z") * rad)
{
}

PHG4GDMLDetector::~PHG4GDMLDetector()
{
}

void

PHG4GDMLDetector::Print(const std::string &what) const
{
  cout << "PHG4GDMLDetector::" << GetName() << " - import " << m_TopVolName << " from " << m_GDMPath << " with shift "
       << m_placeX << ","
       << m_placeY << ","
       << m_placeZ << "cm and rotation "
       << m_rotationX << ","
       << m_rotationY << ","
       << m_rotationZ << "rad" << endl;
}

void PHG4GDMLDetector::Construct(G4LogicalVolume *logicWorld)
{
}
