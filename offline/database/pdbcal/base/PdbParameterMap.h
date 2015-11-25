#ifndef PDBPARAMETERMAP__H
#define PDBPARAMETERMAP__H

#include "PdbCalChan.h"

#include <map>
#include <string>

class PdbParameterMap: public PdbCalChan
{
 public:
  PdbParameterMap() {}
  virtual ~PdbParameterMap() {}

  void print() const;
  void Reset(); // from PHObject - clear content

  std::pair<std::map<const std::string, double>::const_iterator,
    std::map<const std::string, double>::const_iterator> get_dparam_iters() const
  {return make_pair(dparams.begin(),dparams.end());}

  std::pair<std::map<const std::string, int>::const_iterator,
    std::map<const std::string, int>::const_iterator> get_iparam_iters() const
  {return make_pair(iparams.begin(),iparams.end());}

  std::pair<std::map<const std::string, std::string>::const_iterator,
    std::map<const std::string, std::string>::const_iterator> get_cparam_iters() const
  {return make_pair(cparams.begin(),cparams.end());}

  void set_int_param(const std::string &name, const int ival);
  void set_double_param(const std::string &name, const double dval);
  void set_string_param(const std::string &name, const std::string &str);

 protected:

  std::map<const std::string, double> dparams;
  std::map<const std::string, int> iparams;
  std::map<const std::string, std::string> cparams;


  ClassDef(PdbParameterMap,1)

}; 

#endif
