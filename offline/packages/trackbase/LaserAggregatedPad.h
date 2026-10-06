/**
 * @file trackbase/LaserAggregatedPad.h
 * @author Ben Kimelman
 * @date February 2024
 * @brief Base class for laser cluster object
 */
#ifndef TRACKBASE_LASERAGGREGATEDPAD_H
#define TRACKBASE_LASERAGGREGATEDPAD_H

#include <phool/PHObject.h>



#include <iostream>
#include <limits>

/**
 * @brief Base class for laser cluster object
 *
 * Virtual base class for TPC laser cluster object 
 */
class LaserAggregatedPad : public PHObject
{
 public:
  //! dtor
  ~LaserAggregatedPad() override = default;
  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override
  {
    os << "LaserAggregatedPad base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }
  
  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;
  
  //! copy content from base class
  virtual void CopyFrom( const LaserAggregatedPad& /*source*/)  {}

  //! copy content from base class
  virtual void CopyFrom( LaserAggregatedPad* /*source*/)  {}

  virtual void addAdc(double /*adc*/) {}
  virtual double getAdc() const { return std::numeric_limits<double>::max(); }

  virtual void addNHits(double /*nhits*/) {}
  virtual double getNHits() const { return std::numeric_limits<double>::max(); }

 protected:
  LaserAggregatedPad() = default;
  ClassDefOverride(LaserAggregatedPad, 1)
};

#endif //TRACKBASE_LASERAGGREGATEPAD_H
