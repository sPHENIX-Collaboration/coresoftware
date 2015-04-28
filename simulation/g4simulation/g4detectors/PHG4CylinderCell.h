#ifndef PHG4CYLINDERCELL_H
#define PHG4CYLINDERCELL_H

#include <phool/PHObject.h>
#include <cmath>
#include <map>

class PHG4CylinderCell : public PHObject
{
 public:

  virtual ~PHG4CylinderCell(){}

  virtual void identify(std::ostream& os = std::cout) const {
    os << "PHG4CylinderCell base class" << std::endl;
  }
  
  virtual std::pair< std::map<unsigned int,float>::const_iterator, std::map<unsigned int,float>::const_iterator > get_g4hits() {
    std::map <unsigned int, float> dummy;
    return std::make_pair(dummy.begin(), dummy.end());
  }
  
  virtual void add_edep(const unsigned int g4hitid, const float edep) {return;}
  virtual void set_cell_id(const unsigned int id) {return;}
  virtual void set_layer(const unsigned int i) {return;}
 
  virtual double get_edep() const {return NAN;}
  virtual unsigned int get_layer() const {return 0xFFFFFFFF;}
  virtual unsigned int get_cell_id() const {return 0xFFFFFFFF;}
  virtual int get_binz() const {return -1;}
  virtual int get_binphi() const {return -1;}
  virtual int get_bineta() const {return -1;}

  virtual void set_phibin(const int i) {return;}
  virtual void set_zbin(const int i) {return;}
  virtual void set_etabin(const int i) {return;}

  virtual void set_sensor_index(const std::string &si) {return;}
  virtual std::string get_sensor_index() const {return "";}
  virtual void set_ladder_phi_index(int i) {return;}
  virtual int get_ladder_phi_index() const {return -9999;}
  virtual void set_ladder_z_index(int i) {return;}
  virtual int get_ladder_z_index() const {return -9999;}
  
 protected:

  PHG4CylinderCell() {}
  ClassDef(PHG4CylinderCell,1)
};

#endif
