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

/*!
 * \brief PHGenIntegralv1
 */
class PHGenIntegralv1 : public PHGenIntegral
{
 public:
  PHGenIntegralv1();
  explicit PHGenIntegralv1(const std::string& description);
  virtual ~PHGenIntegralv1();

  virtual PHObject* clone() const;
  virtual int isValid() const { return 1; }
  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset();

  virtual int Integrate() const { return 1; }
  /// For integral objects, e.g. integrated luminosity counter, integrate with another object from another run
  virtual int Integrate(PHObject*);
  virtual void CopyContent(PHObject* obj);

  //! Integrated luminosity in pb^-1
  Double_t get_Integrated_Lumi() const
  {
    return fIntegratedLumi;
  }

  //! Integrated luminosity in pb^-1
  void set_Integrated_Lumi(Double_t integratedLumi)
  {
    fIntegratedLumi = integratedLumi;
  }

  //! Number of accepted events in the event generator. This can be higher than fNProcessedEvent depending on trigger on the event generator
  ULong64_t get_N_Generator_Accepted_Event() const
  {
    return fNGeneratorAcceptedEvent;
  }

  //! Number of accepted events in the event generator. This can be higher than fNProcessedEvent depending on trigger on the event generator
  void set_N_Generator_Accepted_Event(ULong64_t nGeneratorAcceptedEvent)
  {
    fNGeneratorAcceptedEvent = nGeneratorAcceptedEvent;
  }

  //! Number of processed events in the Fun4All cycles
  ULong64_t get_N_Processed_Event() const
  {
    return fNProcessedEvent;
  }

  //! Number of processed events in the Fun4All cycles
  void set_N_Processed_Event(ULong64_t nProcessedEvent)
  {
    fNProcessedEvent = nProcessedEvent;
  }

  //! Sum of weight assigned to the events by the generators.
  //! Event weight is normally 1 and thus equal to number of the generated event and is uninteresting.
  //! However, there are several cases where one may have nontrivial event weights, e.g. using user hooks in Pythia8 generators to reweight the phase space
  Double_t get_Sum_Of_Weight() const
  {
    return fSumOfWeight;
  }

  //! Sum of weight assigned to the events by the generators.
  //! Event weight is normally 1 and thus equal to number of the generated event and is uninteresting.
  //! However, there are several cases where one may have nontrivial event weights, e.g. using user hooks in Pythia8 generators to reweight the phase space
  void set_Sum_Of_Weight(Double_t sumOfWeight)
  {
    fSumOfWeight = sumOfWeight;
  }

  //! description on the source
  const std::string& get_Description() const
  {
    return fDescription;
  }

  //! description on the source
  void set_Description(const std::string& description)
  {
    fDescription = description;
  }

 private:
  //! Number of processed events in the Fun4All cycles
  ULong64_t fNProcessedEvent;

  //! Number of accepted events in the event generator. This can be higher than fNProcessedEvent depending on trigger on the event generator
  ULong64_t fNGeneratorAcceptedEvent;

  //! Integrated luminosity in pb^-1
  Double_t fIntegratedLumi;

  //! Sum of weight assigned to the events by the generators.
  //! Event weight is normally 1 and thus equal to number of the generated event and is uninteresting.
  //! However, there are several cases where one may have nontrivial event weights, e.g. using user hooks in Pythia8 generators to reweight the phase space
  Double_t fSumOfWeight;

  //! description on the source
  std::string fDescription;

  ClassDef(PHGenIntegralv1, 1)
};

#endif /* PHHEPMC_PHGENINTEGRALV1_H */
