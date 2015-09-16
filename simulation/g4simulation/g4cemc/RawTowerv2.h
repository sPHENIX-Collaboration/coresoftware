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

  double get_edep() const { return edep; }
  double get_thetaMin() const { return thetaMin; }
  double get_thetaSize() const { return thetaSize; }
  double get_phiMin() const { return phiMin; }
  double get_phiSize() const { return phiSize; }
  double get_zMin() const { return zMin; }
  double get_zSize() const { return zSize; }

  void set_edep( const double set ) { edep = set; }
  void set_thetaMin( const double set ) { thetaMin = set; }
  void set_thetaSize( const double set ) { thetaSize = set; }
  void set_phiMin( const double set ) { phiMin = set; }
  void set_phiSize( const double set ) { phiSize = set; }
  void set_zMin( const double set ) { zMin = set; }
  void set_zSize( const double set ) { zSize = set; }

  RawTowerv2* clone() const;

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
