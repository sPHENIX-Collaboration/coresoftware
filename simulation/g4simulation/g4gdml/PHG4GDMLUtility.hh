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
class PHG4GDMLConfig;
class PHCompositeNode;

/*!
 * \brief PHG4GDMLUtility is utility class that drive the PHG4GDMLWriteStructure
 */
class PHG4GDMLUtility
{
 public:
  virtual ~PHG4GDMLUtility() {}

  //! save the current Geant4 geometry to GDML file. Reading PHG4GDMLConfig from topNode
  static void Dump_GDML(const std::string &filename, G4VPhysicalVolume *vol, PHCompositeNode *topNode = nullptr);

  //! same as above but use default Geant functions as much as possible
  static void Dump_G4_GDML(const std::string &filename, G4VPhysicalVolume *vol);

  static constexpr const char *get_PHG4GDML_Schema()
  {
    return "http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd";
  }

  //! Get or make the parameter node PHG4GDMLConfig from DST nodes. If not found, make a new one
  static PHG4GDMLConfig *GetOrMakeConfigNode(PHCompositeNode *topNode, bool build_new = true);

  static constexpr const char *getDSTNodeName()
  {
    return "G4GDML_CONFIG";
  }

 private:
  PHG4GDMLUtility() {}
};

#endif /* SIMULATION_CORESOFTWARE_SIMULATION_G4SIMULATION_G4GDML_PHG4GDMLUTILITY_HH_ */
