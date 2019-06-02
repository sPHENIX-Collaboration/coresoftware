// $Id: $

/*!
 * \file PHGenIntegralv1.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHGenIntegralv1.h"

#include <phool/PHObject.h>  // for PHObject

#include <cstdlib>
#include <iostream>
using namespace std;

PHGenIntegralv1::PHGenIntegralv1()
{
  Reset();
}

PHGenIntegralv1::PHGenIntegralv1(const std::string& description)
{
  Reset();
  fDescription = description;
}

PHGenIntegralv1::~PHGenIntegralv1()
{
}

PHObject* PHGenIntegralv1::clone() const
{
  //this class is simple, use the default copy constructor for cloning
  return new PHGenIntegralv1(*this);
}

void PHGenIntegralv1::identify(ostream& os) const
{
  os << "PHGenIntegralv1::identify: " << get_Description() << endl
     << " N_Generator_Accepted_Event = " << get_N_Generator_Accepted_Event() << " @ " << get_CrossSection_Generator_Accepted_Event() << " pb" << endl
     << "          N_Processed_Event = " << get_N_Processed_Event() << " @ " << get_CrossSection_Processed_Event() << " pb" << endl
     << "              Sum_Of_Weight = " << get_Sum_Of_Weight() << endl
     << "            Integrated_Lumi = " << get_Integrated_Lumi() << " pb^-1" << endl;
}

void PHGenIntegralv1::Reset()
{
  fNProcessedEvent = 0;
  fNGeneratorAcceptedEvent = 0;
  fIntegratedLumi = 0;
  fSumOfWeight = 0;
  fDescription = "Source Not Provided";
}

int PHGenIntegralv1::Integrate(PHObject* incoming_object)
{
  const PHGenIntegral* in_gen = dynamic_cast<const PHGenIntegral*>(incoming_object);

  if (!in_gen)
  {
    cout << "PHGenIntegralv1::Integrate - Fatal Error - "
         << "input object is not a PHGenIntegral: ";
    incoming_object->identify();

    exit(EXIT_FAILURE);
  }

  if (fIntegratedLumi == 0 and fNProcessedEvent == 0)
  {
    fDescription = in_gen->get_Description();
  }
  else if (fDescription != in_gen->get_Description())
  {
    fDescription = fDescription + ", and " + in_gen->get_Description();
  }

  fNProcessedEvent += in_gen->get_N_Processed_Event();
  fNGeneratorAcceptedEvent += in_gen->get_N_Generator_Accepted_Event();
  fIntegratedLumi += in_gen->get_Integrated_Lumi();
  fSumOfWeight += in_gen->get_Sum_Of_Weight();

  return fNProcessedEvent;
}

void PHGenIntegralv1::CopyContent(PHObject* incoming_object)
{
  const PHGenIntegral* in_gen = dynamic_cast<const PHGenIntegral*>(incoming_object);

  if (!in_gen)
  {
    cout << "PHGenIntegralv1::CopyContent - Fatal Error - "
         << "input object is not a PHGenIntegral: ";
    incoming_object->identify();

    exit(EXIT_FAILURE);
  }

  fNProcessedEvent = in_gen->get_N_Processed_Event();
  fNGeneratorAcceptedEvent = in_gen->get_N_Generator_Accepted_Event();
  fIntegratedLumi = in_gen->get_Integrated_Lumi();
  fSumOfWeight = in_gen->get_Sum_Of_Weight();
  fDescription = in_gen->get_Description();
}
