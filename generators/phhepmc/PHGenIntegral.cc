// $Id: $

/*!
 * \file PHGenIntegral.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHGenIntegral.h"

#include <limits>

//! cross sections for the processed events in pb
Double_t PHGenIntegral::get_CrossSection_Processed_Event() const
{
  return get_Integrated_Lumi() > 0 ? get_N_Processed_Event() / get_Integrated_Lumi() : std::numeric_limits<Double_t>::signaling_NaN();
}

//! cross sections for the events accepted by the event generator in pb
Double_t PHGenIntegral::get_CrossSection_Generator_Accepted_Event() const
{
  return get_Integrated_Lumi() > 0 ? get_N_Generator_Accepted_Event() / get_Integrated_Lumi() : std::numeric_limits<Double_t>::signaling_NaN();
}

//! description on the source
const std::string& PHGenIntegral::get_Description() const
{
  static const std::string s_invalid("Invalid");
  return s_invalid;
}
