// $Id: $

/*!
 * \file PHG4GDMLConfig.hh
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef SIMULATION_CORESOFTWARE_SIMULATION_G4SIMULATION_G4GDML_PHG4GDMLCONFIG_HH_
#define SIMULATION_CORESOFTWARE_SIMULATION_G4SIMULATION_G4GDML_PHG4GDMLCONFIG_HH_

#include <phool/PHObject.h>

#include <iostream>  // for operator<<, basic_ostream, basic_ostream...
#include <set>

class G4VPhysicalVolume;
class G4LogicalVolume;

/*!
 * \brief PHG4GDMLConfig
 */
class PHG4GDMLConfig : public PHObject
{
 public:
  PHG4GDMLConfig() {}
  virtual ~PHG4GDMLConfig() {}

  virtual void Reset()
  {
    excluded_physical_vol.clear();
    excluded_logical_vol.clear();
  }
  virtual int isValid() const { return 1; }
  virtual void identify(std::ostream &os = std::cout) const
  {
    os << "PHG4GDMLConfig with " << excluded_physical_vol.size() << "excluded physical volume and "
              << excluded_logical_vol.size() << " excluded logical volume" << std::endl;
  }

  void exclude_physical_vol(const G4VPhysicalVolume *vol) { excluded_physical_vol.insert(vol); }
  void exclude_logical_vol(const G4LogicalVolume *vol) { excluded_logical_vol.insert(vol); }
  const std::set<const G4VPhysicalVolume *> &get_excluded_physical_vol() const { return excluded_physical_vol; }
  const std::set<const G4LogicalVolume *> &get_excluded_logical_vol() const { return excluded_logical_vol; }

 private:
  std::set<const G4VPhysicalVolume *> excluded_physical_vol;
  std::set<const G4LogicalVolume *> excluded_logical_vol;
};

#endif /* SIMULATION_CORESOFTWARE_SIMULATION_G4SIMULATION_G4GDML_PHG4GDMLCONFIG_HH_ */
