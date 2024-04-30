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
#include <array>

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

  virtual float getLayer() const { return NAN; }
  virtual void setLayer(float) {}
  virtual float getIPhi() const { return NAN; }
  virtual void setIPhi(float) {}
  virtual float getIT() const { return NAN; }
  virtual void setIT(float) {}

  //
  // cluster info
  //
  virtual void setAdc(unsigned int) {}
  virtual unsigned int getAdc() const { return UINT_MAX; }


  virtual unsigned int getNhits() const {return UINT_MAX;}
  virtual void setNhits( unsigned int) {}

  virtual void addHit() {}
 
  virtual void setHitLayer(int, float) {}
  virtual float getHitLayer(int) const { return NAN; }

  virtual void setHitIPhi(int, float) {}
  virtual float getHitIPhi(int) const { return NAN; }

  virtual void setHitIT(int, float) {}
  virtual float getHitIT(int) const { return NAN; }

  virtual void setHitX(int, float) {}
  virtual float getHitX(int) const { return NAN; }

  virtual void setHitY(int, float) {}
  virtual float getHitY(int) const { return NAN; }

  virtual void setHitZ(int, float) {}
  virtual float getHitZ(int) const { return NAN; }
 
  virtual void setHitAdc(int, float) {}
  virtual float getHitAdc(int) const { return NAN; }


 protected:
  LaserCluster() = default;
  ClassDefOverride(LaserCluster, 1)
};

#endif //TRACKBASE_LASERCLUSTER_H
