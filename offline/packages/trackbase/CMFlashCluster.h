/**
 * @file trackbase/CMFlashCluster.h
 * @author Tony Frawley
 * @date January 2022
 * @brief Base class for central membrane flash cluster cluster object
 */
#ifndef TRACKBASE_CMFLASHCLUSTER_H
#define TRACKBASE_CMFLASHCLUSTER_H

#include <phool/PHObject.h>

#include <climits>
#include <cmath>
#include <iostream>
#include <memory>


/**
 * @brief Base class for Central Membrane flassh combined cluster object
 *
 * Virtual base class for TPC CM flash cluster object 
 */
class CMFlashCluster : public PHObject
{
 public:
  //! dtor
  ~CMFlashCluster() override {}
  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override
  {
    os << "CMFlashCluster base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  
  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;
  
  //! copy content from base class
  virtual void CopyFrom( const CMFlashCluster& )  {}

  //! copy content from base class
  virtual void CopyFrom( CMFlashCluster* )  {}
 
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
  virtual unsigned int getNclusters() const {return UINT_MAX;}
  virtual void setNclusters( unsigned int) {}

 protected:
  CMFlashCluster() = default;
  ClassDefOverride(CMFlashCluster, 1)
};

#endif //TRACKBASE_CMFLASHCLUSTER_H
