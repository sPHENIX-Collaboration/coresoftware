/**
 * @file trackbase/LaserCluster.h
 * @author Ben Kimelman
 * @date February 2024
 * @brief Base class for laser cluster object
 */
#ifndef TRACKBASE_LASERCLUSTER_H
#define TRACKBASE_LASERCLUSTER_H

#include <phool/PHObject.h>

#include <climits>
#include <cmath>
#include <iostream>
#include <memory>


/**
 * @brief Base class for laser cluster object
 *
 * Virtual base class for TPC laser cluster object 
 */
class LaserCluster : public PHObject
{
 public:
  //! dtor
  ~LaserCluster() override {}
  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override
  {
    os << "LaserCluster base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  
  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;
  
  //! copy content from base class
  virtual void CopyFrom( const LaserCluster& )  {}

  //! copy content from base class
  virtual void CopyFrom( LaserCluster* )  {}
 
  //
  // cluster position
  //
  virtual float getX() const { return NAN; }
  virtual void setX(float) {}
  virtual float getY() const { return NAN; }
  virtual void setY(float) {}
  virtual float getZ() const { return NAN; }
  virtual void setZ(float) {}

  //
  // cluster info
  //
  virtual void setAdc(unsigned int) {}
  virtual unsigned int getAdc() const { return UINT_MAX; }


  virtual unsigned int getNhits() const {return UINT_MAX;}
  virtual void setNhits( unsigned int) {}


 protected:
  LaserCluster() = default;
  ClassDefOverride(LaserCluster, 1)
};

#endif //TRACKBASE_LASERCLUSTER_H
