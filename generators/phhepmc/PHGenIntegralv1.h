// $Id: $

/*!
 * \file PHGenIntegralv1.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHHEPMC_PHGENINTEGRALV1_H
#define PHHEPMC_PHGENINTEGRALV1_H

#include "PHGenIntegral.h"

#include <iostream>         // for cout, ostream
#include <string>           // for string

class PHObject;

/*!
 * \brief PHGenIntegralv1
 */
class PHGenIntegralv1 : public PHGenIntegral
{
 public:
  PHGenIntegralv1();
  explicit PHGenIntegralv1(const std::string& description);
  ~PHGenIntegralv1() override{}

  PHObject* CloneMe() const override {return new PHGenIntegralv1(*this);}
  int isValid() const override { return 1; }
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;

  int Integrate() const override { return 1; }
  /// For integral objects, e.g. integrated luminosity counter, integrate with another object from another run
  int Integrate(PHObject*) override;
  void CopyFrom(const PHObject* obj) override;

  //! Integrated luminosity in pb^-1
  Double_t get_Integrated_Lumi() const override
  {
    return m_IntegratedLumi;
  }

  //! Integrated luminosity in pb^-1
  void set_Integrated_Lumi(Double_t integratedLumi) override
  {
    m_IntegratedLumi = integratedLumi;
  }

  //! Number of accepted events in the event generator. This can be higher than m_NProcessedEvent depending on trigger on the event generator
  ULong64_t get_N_Generator_Accepted_Event() const override
  {
    return m_NGeneratorAcceptedEvent;
  }

  //! Number of accepted events in the event generator. This can be higher than m_NProcessedEvent depending on trigger on the event generator
  void set_N_Generator_Accepted_Event(ULong64_t nGeneratorAcceptedEvent) override
  {
    m_NGeneratorAcceptedEvent = nGeneratorAcceptedEvent;
  }

  //! Number of processed events in the Fun4All cycles
  ULong64_t get_N_Processed_Event() const override
  {
    return m_NProcessedEvent;
  }

  //! Number of processed events in the Fun4All cycles
  void set_N_Processed_Event(ULong64_t nProcessedEvent) override
  {
    m_NProcessedEvent = nProcessedEvent;
  }

  //! Sum of weight assigned to the events by the generators.
  //! Event weight is normally 1 and thus equal to number of the generated event and is uninteresting.
  //! However, there are several cases where one may have nontrivial event weights, e.g. using user hooks in Pythia8 generators to reweight the phase space
  Double_t get_Sum_Of_Weight() const override
  {
    return m_SumOfWeight;
  }

  //! Sum of weight assigned to the events by the generators.
  //! Event weight is normally 1 and thus equal to number of the generated event and is uninteresting.
  //! However, there are several cases where one may have nontrivial event weights, e.g. using user hooks in Pythia8 generators to reweight the phase space
  void set_Sum_Of_Weight(Double_t sumOfWeight) override
  {
    m_SumOfWeight = sumOfWeight;
  }

  //! description on the source
  const std::string& get_Description() const override
  {
    return m_Description;
  }

  //! description on the source
  void set_Description(const std::string& description) override
  {
    m_Description = description;
  }

 private:
  //! Number of processed events in the Fun4All cycles
  ULong64_t m_NProcessedEvent;

  //! Number of accepted events in the event generator. This can be higher than m_NProcessedEvent depending on trigger on the event generator
  ULong64_t m_NGeneratorAcceptedEvent;

  //! Integrated luminosity in pb^-1
  Double_t m_IntegratedLumi;

  //! Sum of weight assigned to the events by the generators.
  //! Event weight is normally 1 and thus equal to number of the generated event and is uninteresting.
  //! However, there are several cases where one may have nontrivial event weights, e.g. using user hooks in Pythia8 generators to reweight the phase space
  Double_t m_SumOfWeight;

  //! description on the source
  std::string m_Description;

  ClassDefOverride(PHGenIntegralv1, 1)
};

#endif /* PHHEPMC_PHGENINTEGRALV1_H */
