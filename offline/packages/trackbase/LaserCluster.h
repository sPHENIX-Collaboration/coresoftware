/**
 * @file trackbase/LaserCluster.h
 * @author Ben Kimelman
 * @date February 2024
 * @brief Base class for laser cluster object
 */
#ifndef TRACKBASE_LASERCLUSTER_H
#define TRACKBASE_LASERCLUSTER_H

#include <phool/PHObject.h>

#include <iostream>
#include <limits>

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
  virtual float getX() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void setX(float) {}
  virtual float getY() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void setY(float) {}
  virtual float getZ() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void setZ(float) {}

  virtual bool getFitMode() const { return false; }
  virtual void setFitMode(bool) {}

  virtual float getLayer() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void setLayer(float) {}
  virtual float getIPhi() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void setIPhi(float) {}
  virtual float getIT() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void setIT(float) {}

  //
  // cluster info
  //
  virtual void setAdc(unsigned int) {}
  virtual unsigned int getAdc() const { return std::numeric_limits<unsigned int>::max(); }


  virtual unsigned int getNhits() const {return std::numeric_limits<unsigned int>::max();}
  virtual void setNhits( unsigned int) {}

  virtual void setNLayers(unsigned int) {}
  virtual unsigned int getNLayers() const { return std::numeric_limits<unsigned int>::max();}

  virtual void setNIPhi(unsigned int) {}
  virtual unsigned int getNIPhi() const { return std::numeric_limits<unsigned int>::max();}

  virtual void setNIT(unsigned int) {}
  virtual unsigned int getNIT() const { return std::numeric_limits<unsigned int>::max();}

  virtual void setSDLayer(float) {}
  virtual float getSDLayer() const { return std::numeric_limits<float>::quiet_NaN(); }

  virtual void setSDIPhi(float) {}
  virtual float getSDIPhi() const { return std::numeric_limits<float>::quiet_NaN(); }

  virtual void setSDIT(float) {}
  virtual float getSDIT() const { return std::numeric_limits<float>::quiet_NaN(); }

  virtual void setSDWeightedLayer(float) {}
  virtual float getSDWeightedLayer() const { return std::numeric_limits<float>::quiet_NaN(); }

  virtual void setSDWeightedIPhi(float) {}
  virtual float getSDWeightedIPhi() const { return std::numeric_limits<float>::quiet_NaN(); }

  virtual void setSDWeightedIT(float) {}
  virtual float getSDWeightedIT() const { return std::numeric_limits<float>::quiet_NaN(); }


  virtual void addHit() {}
 
  virtual void setHitLayer(int, float) {}
  virtual float getHitLayer(int) const { return std::numeric_limits<float>::quiet_NaN(); }

  virtual void setHitIPhi(int, float) {}
  virtual float getHitIPhi(int) const { return std::numeric_limits<float>::quiet_NaN(); }

  virtual void setHitIT(int, float) {}
  virtual float getHitIT(int) const { return std::numeric_limits<float>::quiet_NaN(); }

  virtual void setHitX(int, float) {}
  virtual float getHitX(int) const { return std::numeric_limits<float>::quiet_NaN(); }

  virtual void setHitY(int, float) {}
  virtual float getHitY(int) const { return std::numeric_limits<float>::quiet_NaN(); }

  virtual void setHitZ(int, float) {}
  virtual float getHitZ(int) const { return std::numeric_limits<float>::quiet_NaN(); }
 
  virtual void setHitAdc(int, float) {}
  virtual float getHitAdc(int) const { return std::numeric_limits<float>::quiet_NaN(); }


 protected:
  LaserCluster() = default;
  ClassDefOverride(LaserCluster, 1)
};

#endif //TRACKBASE_LASERCLUSTER_H
