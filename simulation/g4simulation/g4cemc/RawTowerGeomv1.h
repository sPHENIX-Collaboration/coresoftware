#ifndef RawTowerGeomv1_H__
#define RawTowerGeomv1_H__

#include "RawTowerGeom.h"

#include <map>
#include <string>

class RawTowerGeomv1: public RawTowerGeom
{
 public:
  RawTowerGeomv1();

  virtual ~RawTowerGeomv1() {}

  void identify(std::ostream& os = std::cout) const;
  double get_radius() const {return radius;}
  double get_thickness() const {return thickness;}
  int get_phibins() const {return  nphibins;}
  double get_phistep() const {return phistep;}
  double get_phimin() const {return phimin;}
  int get_etabins() const {return netabins;}
  double get_etastep() const {return etastep;}
  double get_etamin() const {return etamin;}

  std::pair<double, double> get_phibounds(const int ibin) const;
  std::pair<double, double> get_etabounds(const int ibin) const;
  double get_etacenter(const int ibin) const;
  double get_phicenter(const int ibin) const;

  int get_etabin(const double eta) const;
  int get_phibin(const double phi) const;

   void set_radius(const double r) {radius = r;}
   void set_thickness(const double t) {thickness = t;}
   void set_phibins(const int i) {nphibins = i;}
   void set_phistep(const double phi) {phistep = phi;}
   void set_phimin(const double phi) {phimin = phi;}
   void set_etabins(const int i) {netabins = i;}
   void set_etamin(const double z) {etamin = z;}
   void set_etastep(const double z) {etastep = z;}
  
 protected:
  double radius;
  double thickness;
  int netabins;
  double etamin;
  double etastep;
  int nphibins;
  double phimin;
  double phistep;
  ClassDef(RawTowerGeomv1,1)
};

#endif
