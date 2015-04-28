#ifndef RAWTOWERV2_H_
#define RAWTOWERV2_H_

#include "RawTowerv1.h"

#include <map>

/**
 * \brief Calorimeter tower including basic geometry information.
 *
 * \author Nils Feege <nils.feege@stonybrook.edu>
 *
 */
class RawTowerv2 : public RawTowerv1 {

public:
  RawTowerv2();
  RawTowerv2(const int ieta, const int iphi);
  virtual ~RawTowerv2();

  void identify(std::ostream& os=std::cout) const;

  int get_bintheta() const { return bineta; }

  double get_edep() { return edep; }
  double get_thetaMin() { return thetaMin; }
  double get_thetaSize() { return thetaSize; }
  double get_phiMin() { return phiMin; }
  double get_phiSize() { return phiSize; }
  double get_zMin() { return zMin; }
  double get_zSize() { return zSize; }

  void set_edep( double set ) { edep = set; }
  void set_thetaMin( double set ) { thetaMin = set; }
  void set_thetaSize( double set ) { thetaSize = set; }
  void set_phiMin( double set ) { phiMin = set; }
  void set_phiSize( double set ) { phiSize = set; }
  void set_zMin( double set ) { zMin = set; }
  void set_zSize( double set ) { zSize = set; }

  RawTowerv2* clone();

protected:

  int bineta;
  int binphi;

  double thetaMin;
  double thetaSize;
  double phiMin;
  double phiSize;
  double zMin;
  double zSize;
  double edep;

  ClassDef(RawTowerv2,1)
};

#endif /* RAWTOWERV1_H_ */
