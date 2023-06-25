// $Id: $

/*!
 * \file RawCluster.cc
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \author Francesco Vassalli <Francesco.Vassalli@colorado.edu>
 * \version $Revision: 2  $
 * \date $Date: 07/13/18$
 */

#include "RawCluster.h"

#include <cstdlib>
#include <iostream>

using namespace std;

std::pair<const std::string, RawCluster::PROPERTY_TYPE>
RawCluster::get_property_info(const PROPERTY prop_id)
{
  switch (prop_id)
  {
  case prop_ecore:
    return make_pair("EM cluster energy estimated with its core towers", RawCluster::type_float);
    break;
  case prop_prob:
    return make_pair("cluster template probability for EM shower", RawCluster::type_float);
    break;
  case prop_chi2:
    return make_pair("reduced chi2 for EM shower", RawCluster::type_float);
    break;

  case prop_et_iso_calotower_sub_R01:
    return make_pair("subtracted calortower isolation ET R=.1", RawCluster::type_float);
    break;

  case prop_et_iso_calotower_sub_R02:
    return make_pair("subtracted calortower isolation ET R=.2", RawCluster::type_float);
    break;
  case prop_et_iso_calotower_sub_R03:
    return make_pair("subtracted calortower isolation ET R=.3", RawCluster::type_float);
    break;
  case prop_et_iso_calotower_sub_R04:
    return make_pair("subtracted calortower isolation ET R=.4", RawCluster::type_float);
    break;
  case prop_et_iso_calotower_R01:
    return make_pair("calortower isolation ET R=.1", RawCluster::type_float);
    break;
  case prop_et_iso_calotower_R02:
    return make_pair("calortower isolation ET R=.2", RawCluster::type_float);

    break;
  case prop_et_iso_calotower_R03:
    return make_pair("calortower isolation ET R=.3", RawCluster::type_float);
    break;
  case prop_et_iso_calotower_R04:
    return make_pair("calortower isolation ET R=.4", RawCluster::type_float);
    break;
    //  case prop_truth_track_ID:
    //    return make_pair("truth cluster's PHG4Particle ID", RawCluster::type_int);
    //    break;
    //  case prop_truth_flavor:
    //    return make_pair("truth cluster's PHG4Particle flavor", RawCluster::type_int);
    //    break;

  default:
    cout << "RawCluster::get_property_info - Fatal Error - unknown index " << prop_id << endl;
    exit(1);
  }
}

bool RawCluster::check_property(const PROPERTY prop_id, const PROPERTY_TYPE prop_type)
{
  pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
  if (property_info.second != prop_type)
  {
    return false;
  }
  return true;
}

string
RawCluster::get_property_type(const PROPERTY_TYPE prop_type)
{
  switch (prop_type)
  {
  case type_int:
    return "int";
  case type_uint:
    return "unsigned int";
  case type_float:
    return "float";
  default:
    return "unkown";
  }
}
