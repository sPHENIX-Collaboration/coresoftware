/**
 * @file trackbase/LaserAggregatedPadv1.h
 * @author Ben Kimelman
 * @date July 2026
 * @brief Version 4 of CMFLashCluster
 */
#ifndef TRACKBASE_LASERAGGREGATEDPADV1_H
#define TRACKBASE_LASERAGGREGATEDPADV1_H

#include "LaserAggregatedPad.h"

#include <iostream>
#include <vector>

class PHObject;

/**
 * @brief Version 1 of LaserAggregatedPad
 *
 * Note - D. McGlinchey June 2018:
 *   CINT does not like "override", so ignore where CINT
 *   complains. Should be checked with ROOT 6 once
 *   migration occurs.
 */


class LaserAggregatedPadv1 : public LaserAggregatedPad
{
 public:
  //! ctor
  LaserAggregatedPadv1() = default;

  // PHObject virtual overloads
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new LaserAggregatedPadv1(*this); }
 
  //! copy content from base class
  void CopyFrom( const LaserAggregatedPad& ) override;

  //! copy content from base class
  void CopyFrom( LaserAggregatedPad* source ) override
  { CopyFrom( *source ); }

  void addAdc(double adc) override { m_adc += adc; }
  double getAdc() const override { return m_adc; }

  void addNHits(double nhits) override { m_nHits += nhits; }
  double getNHits() const override { return m_nHits; }

  void identify(std::ostream& os = std::cout) const override;

 private:

  double m_adc {0.0};
  double m_nHits {0.0};

  ClassDefOverride(LaserAggregatedPadv1, 1)
};

#endif //TRACKBASE_LASERAGGREGATEDPADV1_H
