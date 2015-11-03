// $Id: PdbCalParameters.hh,v 1.1 2014/12/31 18:56:10 jinhuang Exp $                                                                                             

/*!
 * \file PdbCalParameters.h
 * \brief Generic purpose calibration parameter sets
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.1 $
 * \date $Date: 2014/12/31 18:56:10 $
 */

#ifndef PDBCALPARAMETERS_H_
#define PDBCALPARAMETERS_H_

#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include "PdbCalChan.hh"

/*!
 * \brief Generic purpose calibration parameter sets
 */
class PdbCalParameters : public PdbCalChan
{
public:
  PdbCalParameters();
  virtual
  ~PdbCalParameters();

  virtual void
  print() const;

  size_t size() const
  {
    return _data.size();
  }

  bool
  is_parameter_exist(const std::string & name) const
  {
    Data_t::const_iterator it = _data.find(name);

    return it != _data.end();
  }

  template<class Type>
    Type
    get_parameter(const std::string & name) const
    {

      Type value = Type();

      Data_t::const_iterator it = _data.find(name);

      if (it == _data.end())
        {
          std::cout
              << "PdbCalParameters::get_value - Error - can not find parameter with name "
              << name << ". Returning " << value << ". Available entries are:"
              << std::endl;
          print();
          return value;
        }

      if (it->second.length() == 0)
        return value;

      std::istringstream i(it->second);

      i >> value;

      if (i.fail())
        {
          std::cout
              << "PdbCalParameters::get_value - Error - can not read parameter with name "
              << name << ". Returning " << value << ". Available entries are:"
              << std::endl;
          print();
        }

      return value;
    }

  template<class Type>
    void
    set_parameter(const std::string & name, const Type & value)
    {
      std::ostringstream o;

      o << value;

      _data[name] = o.str();
    }

private:

  typedef std::map<std::string, std::string> Data_t;

  Data_t _data;

ClassDef(PdbCalParameters,1)
};

#endif /* PDBCALPARAMETERS_H_ */
