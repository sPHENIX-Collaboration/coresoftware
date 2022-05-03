// $Id: $

/*!
 * \file PHGenIntegral.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHHEPMC_PHGENINTEGRAL_H
#define PHHEPMC_PHGENINTEGRAL_H

#include <phool/PHObject.h>

#include <string>

/*!
 * \brief PHGenIntegral
 */
class PHGenIntegral : public PHObject
{
 public:
  PHGenIntegral(){}
  ~PHGenIntegral() override{}

  //! Integrated luminosity in pb^-1
  virtual Double_t get_Integrated_Lumi() const
  {
    return 0;
  }

  //! Integrated luminosity in pb^-1
  virtual void set_Integrated_Lumi(Double_t /*integratedLumi*/)
  {
  }

  //! Number of accepted events in the event generator. This can be higher than fNProcessedEvent depending on trigger on the event generator
  virtual ULong64_t get_N_Generator_Accepted_Event() const
  {
    return 0;
  }

  //! Number of accepted events in the event generator. This can be higher than fNProcessedEvent depending on trigger on the event generator
  virtual void set_N_Generator_Accepted_Event(ULong64_t /*nGeneratorAcceptedEvent*/)
  {
  }

  //! Number of processed events in the Fun4All cycles
  virtual ULong64_t get_N_Processed_Event() const
  {
    return 0;
  }

  //! Number of processed events in the Fun4All cycles
  virtual void set_N_Processed_Event(ULong64_t /*nProcessedEvent*/)
  {
  }

  //! Sum of weight assigned to the events by the generators.
  //! Event weight is normally 1 and thus equal to number of the generated event and is uninteresting.
  //! However, there are several cases where one may have nontrivial event weights, e.g. using user hooks in Pythia8 generators to reweight the phase space
  virtual Double_t get_Sum_Of_Weight() const
  {
    return 0;
  }

  //! Sum of weight assigned to the events by the generators.
  //! Event weight is normally 1 and thus equal to number of the generated event and is uninteresting.
  //! However, there are several cases where one may have nontrivial event weights, e.g. using user hooks in Pythia8 generators to reweight the phase space
  virtual void set_Sum_Of_Weight(Double_t /*sumOfWeight*/)
  {
  }

  //! cross sections for the processed events in pb
  virtual Double_t get_CrossSection_Processed_Event() const;

  //! cross sections for the events accepted by the event generator in pb
  virtual Double_t get_CrossSection_Generator_Accepted_Event() const;

  //! description on the source
  virtual const std::string& get_Description() const;

  //! description on the source
  virtual void set_Description(const std::string& /*description*/)
  {
  }

  ClassDefOverride(PHGenIntegral, 1)
};

#endif /* PHHEPMC_PHGENINTEGRAL_H */
