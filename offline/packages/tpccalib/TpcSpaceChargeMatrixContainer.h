#ifndef TPCCALIB_TPCSPACECHARGEMATRIXCONTAINER_H
#define TPCCALIB_TPCSPACECHARGEMATRIXCONTAINER_H

/**
 * @file tpccalib/TpcSpaceChargeMatrixContainer.h
 * @author Hugo Pereira Da Costa
 * @date June 2018
 * @brief Contains matrices needed for space charge trackbase reconstruction
 */

#include <phool/PHObject.h>

#include <array>

/**
 * @brief Cluster container object
 */
class TpcSpaceChargeMatrixContainer : public PHObject
{
  public:

  /// constructor
  TpcSpaceChargeMatrixContainer()
  {}

  /// destructor
  ~TpcSpaceChargeMatrixContainer() override = default;

  ///@name accessors
  //@{

  /// identify object
  void identify(std::ostream &/*os*/ = std::cout) const override
  {}

  /// get grid dimensions
  virtual void get_grid_dimensions( int& /*phibins*/, int& /*rbins*/, int& /*zbins*/ ) const
  {}

  /// get total grid size
  virtual int get_grid_size() const
  { return 0; }

  /// get grid index for given sub-indexes
  virtual int get_cell_index( int /*iphibin*/, int /*irbin*/, int /*izbin*/ ) const
  { return -1; }

  /// get entries for a given cell
  virtual int get_entries( int /*cell_index*/ ) const
  { return 0; }

  /// get left hand side
  virtual float get_lhs( int /*cell_index*/, int /*i*/, int /*j*/ ) const
  { return 0; }

  /// get right hand side
  virtual float get_rhs( int /*cell_index*/, int /*i*/ ) const
  { return 0; }

  //@}

  ///@name modifiers
  //@{

  /// reset method
  void Reset() override
  {}

  /// set grid dimensions
  /**
  \param phibins the number of bins in the azimuth direction
  \param zbins the number of bins along z
  */
  virtual void set_grid_dimensions( int /*phibins*/, int /*rbins*/, int /*zbins*/ )
  {}

  /// increment cell entries
  virtual void add_to_entries( int /*cell_index*/ )
  {}

  /// increment cell entries
  virtual void add_to_entries( int /*cell_index*/, int /*value*/ )
  {}

  /// increment left hand side matrix
  virtual void add_to_lhs( int /*cell_index*/, int /*i*/, int /*j*/, float /*value*/ )
  {}

  /// increment right hand side column
  virtual void add_to_rhs( int /*cell_index*/, int /*i*/, float /*value*/ )
  {}

  /// add content from other container, returns true on success
  virtual bool add( const TpcSpaceChargeMatrixContainer& /*other*/ )
  { return false; }

  //@}

  private:

  ClassDefOverride(TpcSpaceChargeMatrixContainer, 1)

};

#endif
