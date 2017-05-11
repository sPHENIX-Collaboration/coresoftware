// $Id: $

/*!
 * \file PHG4GDMLUtility.hh
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef CORESOFTWARE_SIMULATION_G4SIMULATION_G4GDML_PHG4GDMLUTILITY_HH_
#define CORESOFTWARE_SIMULATION_G4SIMULATION_G4GDML_PHG4GDMLUTILITY_HH_

#include <string>

class G4VPhysicalVolume;

/*!
 * \brief PHG4GDMLUtility is utility class that drive the PHG4GDMLWriteStructure
 */
class PHG4GDMLUtility
{
 public:
  virtual ~PHG4GDMLUtility();

  static void Dump_GDML(const std::string &filename, G4VPhysicalVolume * vol);

  static constexpr const char * get_PHG4GDML_Schema()
  {
    return "http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd";
  }

 private:
  PHG4GDMLUtility();
};

#endif /* SIMULATION_CORESOFTWARE_SIMULATION_G4SIMULATION_G4GDML_PHG4GDMLUTILITY_HH_ */
