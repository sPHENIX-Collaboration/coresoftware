
#ifndef PHBFIELDMAP_HH
#define PHBFIELDMAP_HH

#include "phool/phool.h"
#include "TObject.h"

class PHBFieldMap : public TObject
{

public:


  PHBFieldMap(const std::string& mapfile);
  virtual ~PHBFieldMap();

  // fetch b field value, interface to outside
  bool get_bfield( const double *loc, double *field );

  double get_scale_factor(){return _scale_factor;}
  void set_scale_factor(double new_factor);

private:

  void *_fieldmap;
  double _scale_factor; 

};

#endif


