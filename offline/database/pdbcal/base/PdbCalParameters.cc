// $Id: PdbCalParameters.cc,v 1.1 2014/12/31 18:56:10 jinhuang Exp $                                                                                             

/*!
 * \file PdbCalParameters.cc
 * \brief Generic purpose calibration parameter sets
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.1 $
 * \date $Date: 2014/12/31 18:56:10 $
 */

#include "PdbCalParameters.hh"

using namespace std;

PdbCalParameters::PdbCalParameters()
{
  // TODO Auto-generated constructor stub

}

PdbCalParameters::~PdbCalParameters()
{
  // TODO Auto-generated destructor stub
}

void
PdbCalParameters::print() const
{
  cout << "PdbCalParameters::print - contains " << _data.size() << " entries:"
      << endl;

  for (Data_t::const_iterator it1 = _data.begin(); it1 != _data.end(); ++it1)
    {
      cout << "|--- " << (*it1).first << " = " << (*it1).second << endl;
    }
}
