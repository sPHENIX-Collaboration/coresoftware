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
  : m_NProcessedEvent(0)
  , m_NGeneratorAcceptedEvent(0)
  , m_IntegratedLumi(0.)
  , m_SumOfWeight(0.)
  , m_Description("Source Not Provided")
{}

PHGenIntegralv1::PHGenIntegralv1(const std::string& description)
  : m_NProcessedEvent(0)
  , m_NGeneratorAcceptedEvent(0)
  , m_IntegratedLumi(0.)
  , m_SumOfWeight(0.)
  , m_Description(description)
{}

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
  m_NProcessedEvent = 0;
  m_NGeneratorAcceptedEvent = 0;
  m_IntegratedLumi = 0;
  m_SumOfWeight = 0;
  m_Description = "Source Not Provided";
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

  if (m_IntegratedLumi == 0 and m_NProcessedEvent == 0)
  {
    m_Description = in_gen->get_Description();
  }
  else if (m_Description != in_gen->get_Description())
  {
    m_Description = m_Description + ", and " + in_gen->get_Description();
  }

  m_NProcessedEvent += in_gen->get_N_Processed_Event();
  m_NGeneratorAcceptedEvent += in_gen->get_N_Generator_Accepted_Event();
  m_IntegratedLumi += in_gen->get_Integrated_Lumi();
  m_SumOfWeight += in_gen->get_Sum_Of_Weight();

  return m_NProcessedEvent;
}

void PHGenIntegralv1::CopyFrom(const PHObject* incoming_object)
{
  const PHGenIntegral* in_gen = dynamic_cast<const PHGenIntegral*>(incoming_object);

  if (!in_gen)
  {
    cout << "PHGenIntegralv1::CopyFrom - Fatal Error - "
         << "input object is not a PHGenIntegral: ";
    incoming_object->identify();

    exit(EXIT_FAILURE);
  }

  m_NProcessedEvent = in_gen->get_N_Processed_Event();
  m_NGeneratorAcceptedEvent = in_gen->get_N_Generator_Accepted_Event();
  m_IntegratedLumi = in_gen->get_Integrated_Lumi();
  m_SumOfWeight = in_gen->get_Sum_Of_Weight();
  m_Description = in_gen->get_Description();
}
