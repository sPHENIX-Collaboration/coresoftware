/**
 * @file outertracker/OuterTrackerHit.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief OuterTracker hit object
 */
#ifndef OTRACK_OUTERTRACKERHIT_H
#define OTRACK_OUTERTRACKERHIT_H

#include <trackbase/TrkrHit.h>

#include <iostream>

/**
 * @brief OuterTracker hit object
 *
 * Container for OuterTracker hit object representing
 * a single hit pixel within a layer (one layer per instance
 */
class OuterTrackerHit : public TrkrHit
{
 public:
  //! ctor
  OuterTrackerHit();
  //! dtor
  virtual ~OuterTrackerHit(){};

  // PHObject virtual overloads
  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset();
  virtual int isValid() const;

  //void addEnergy(const double edep) {m_edep += edep; std::cout << "added energy " << edep << " to OuterTrackerHit" << std::endl;}
  //double getEnergy() const {return m_edep;}


  /**
   * @brief Set the ADC information
   * @param[in] adc ADC value
   */
  //void setAdc(const short adc) { m_adc = adc; }

  /**
   * @brief Get ADC value
   * @param[out] ADC value
   */
  //short getAdc() const { return m_adc; }

 private:
   //short m_adc;
   //double m_edep;
  ClassDef(OuterTrackerHit, 1);
};

#endif  //OTRACK_OUTERTRACKERHIT_H
