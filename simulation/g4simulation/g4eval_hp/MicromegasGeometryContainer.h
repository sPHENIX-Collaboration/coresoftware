#ifndef G4EVAL_HP_MICROMEGASGEOMETRYCONTAINER_H
#define G4EVAL_HP_MICROMEGASGEOMETRYCONTAINER_H

/**
 * @file tpccalib/MicromegasGeometryContainer.h
 * @author Hugo Pereira Da Costa
 * @date June 2018
 * @brief Contains micromegas strip positions
 */

#include <phool/PHObject.h>
#include <TVector3.h>
#include <map>

/**
 * @brief Cluster container object
 */
class MicromegasGeometryContainer : public PHObject
{
  public:

  /// constructor
  MicromegasGeometryContainer() = default;

  /// destructor
  ~MicromegasGeometryContainer() override = default;

  ///@name accessors
  //@{

  /// identify object
  void identify(std::ostream &/*os*/ = std::cout) const override;

  /// get strip begin from layer, tile and strip number
  TVector3 get_strip_begin( unsigned int /*layer*/, unsigned int /*tile*/, unsigned int /*strip*/ ) const;

  /// get strip begin from layer, tile and strip number
  TVector3 get_strip_end( unsigned int /*layer*/, unsigned int /*tile*/, unsigned int /*strip*/ ) const;
  
  //@}

  ///@name modifiers
  //@{
  
  /// reset method
  void Reset() override;

  /// add strip information
  void add_strip( unsigned int /*layer*/, unsigned int /*tile*/, unsigned int /*strip*/, const TVector3& /*begin*/, const TVector3& /*end*/ );

  //@}
  
  class StripId
  {
    public:
    unsigned int layer = 0;
    unsigned int tile = 0;
    unsigned int strip = 0;
    
    bool operator == (const StripId& other ) const
    { return other.layer == layer && other.tile == tile && other.strip == strip; }
    
    bool operator < (const StripId& other ) const
    {
      if( layer != other.layer ) return layer < other.layer;
      else if( tile != other.tile ) return tile < other.tile;
      else return strip < other.strip;
    }
  };

  private:

  using strip_map_t = std::map<StripId, TVector3>;
  strip_map_t m_strip_begin;
  strip_map_t m_strip_end;
  
  ClassDefOverride(MicromegasGeometryContainer, 1)

};

#endif
