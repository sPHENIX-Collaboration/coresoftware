#ifndef RAWTOWER_PROTOTYPE3_H_
#define RAWTOWER_PROTOTYPE3_H_

#include <calobase/RawTower.h>
#include <calobase/RawTowerDefs.h>
#include <map>
#include <stdint.h>

#include "PROTOTYPE3_FEM.h"

class RawTower_Prototype3 : public RawTower {
 public:
  RawTower_Prototype3();
  RawTower_Prototype3(const RawTower& tower);
  RawTower_Prototype3(RawTowerDefs::keytype id);
  RawTower_Prototype3(const unsigned int icol, const unsigned int irow);
  RawTower_Prototype3(const RawTowerDefs::CalorimeterId caloid, const unsigned int ieta,
             const unsigned int iphi);
  virtual ~RawTower_Prototype3();

  void Reset();
  int isValid() const;
  void identify(std::ostream& os = std::cout) const;

  void set_id(RawTowerDefs::keytype id) { towerid = id; }
  RawTowerDefs::keytype get_id() const { return towerid; }
  int get_bineta() const { return RawTowerDefs::decode_index1(towerid); }
  int get_binphi() const { return RawTowerDefs::decode_index2(towerid); }
  int get_column() const { return RawTowerDefs::decode_index1(towerid); }
  int get_row() const { return RawTowerDefs::decode_index2(towerid); }
  double get_energy() const { return energy; }
  void set_energy(const double e) { energy = e; }
  float get_time() const { return time; }
  void set_time(const float t) { time = t; }

  //---Raw data access------------------------------------------------------------

  enum
  {
    NSAMPLES = PROTOTYPE3_FEM::NSAMPLES
  };
  typedef float signal_type;

  void set_signal_samples(int i,signal_type sig);
  signal_type get_signal_samples(int i) const;
  void set_HBD_channel_number(int i)
    { HBD_channel=i; }
  int get_HBD_channel_number() const
    { return HBD_channel; }

  //---Fits------------------------------------------------------------

    double get_energy_power_law_exp(int verbosity = 0);


 protected:
  RawTowerDefs::keytype towerid;

  //! energy assigned to the tower. Depending on stage of process and DST node
  //! name, it could be energy deposition, light yield or calibrated energies
  double energy;
  //! Time stamp assigned to the tower. Depending on the tower maker, it could
  //! be rise time or peak time.
  float time;

  //Signal samples from DATA
  signal_type signal_samples[NSAMPLES];  //Low Gain
  int HBD_channel;

  ClassDef(RawTower_Prototype3, 3)
};

#endif /* RAWTOWER_PROTOTYPE3_H_ */
