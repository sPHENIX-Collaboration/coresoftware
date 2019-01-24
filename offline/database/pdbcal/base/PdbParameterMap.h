#ifndef PDBCAL_BASE_PDBPARAMETERMAP_H
#define PDBCAL_BASE_PDBPARAMETERMAP_H

#include "PdbCalChan.h"

#include <map>
#include <string>

class PdbParameterMap: public PdbCalChan
{
 public:
  typedef std::map<const std::string, double> dMap;
  typedef std::map<const std::string, int> iMap;
  typedef std::map<const std::string, std::string> strMap;
  typedef dMap::const_iterator dIter;
  typedef iMap::const_iterator iIter;
  typedef strMap::const_iterator strIter;
  typedef std::pair<dIter, dIter> dConstRange;
  typedef std::pair<iIter, iIter> iConstRange;
  typedef std::pair<strIter, strIter> strConstRange;

  PdbParameterMap() {}
  virtual ~PdbParameterMap() {}

  void print() const;
  void Reset(); // from PHObject - clear content

  //! hash of binary information for checking purpose
  size_t get_hash() const;

  dConstRange get_dparam_iters() const 
    {return make_pair(dparams.begin(),dparams.end());}

  iConstRange get_iparam_iters() const 
    {return make_pair(iparams.begin(),iparams.end());}


  strConstRange get_cparam_iters() const 
  {return make_pair(cparams.begin(),cparams.end());}

  void set_int_param(const std::string &name, const int ival);
  void set_double_param(const std::string &name, const double dval);
  void set_string_param(const std::string &name, const std::string &str);

 protected:

  dMap dparams;
  iMap iparams;
  strMap cparams;


  ClassDef(PdbParameterMap,1)

}; 

#endif
