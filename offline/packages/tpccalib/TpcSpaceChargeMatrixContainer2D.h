#ifndef TPCCALIB_TPCSPACECHARGEMATRIXCONTAINERV3_H
#define TPCCALIB_TPCSPACECHARGEMATRIXCONTAINERV3_H

/**
 * @file tpccalib/TpcSpaceChargeMatrixContainer.h
 * @author Xudong Yu
 * @date Februray 2025
 * @brief Contains matrices needed for space charge trackbase reconstruction 2D
 */

#include "TpcSpaceChargeMatrixContainer.h"

#include <array>

/**
 * @brief Cluster container object
 */
class TpcSpaceChargeMatrixContainer2D : public TpcSpaceChargeMatrixContainer
{
 public:
  /// constructor
  TpcSpaceChargeMatrixContainer2D();

  /// destructor
  ~TpcSpaceChargeMatrixContainer2D() override = default;

  ///@name accessors
  //@{

  /// identify object
  void identify(std::ostream& os = std::cout) const override;

  /// get grid dimensions
  void get_grid_dimensions(int& pbins, int& rbins, int& zbins) const override;

  /// get grid size
  int get_grid_size() const override;

  /// get grid index for given sub-indexes
  int get_cell_index(int irbin, int izbin) const override;

  /// get entries for a given cell
  int get_entries(int cell_index) const override;

  /// get left hand side
  float get_lhs(int cell_index, int i, int j) const override;

  /// get right hand side
  float get_rhs(int cell_index, int i) const override;

  //@}

  ///@name modifiers
  //@{

  /// reset method
  void Reset() override;

  /// set grid dimensions
  /**
  \param rbins the number of bins along r
  \param zbins the number of bins along z
  */
  void set_grid_dimensions(int pbins, int rbins, int zbins) override;

  /// increment cell entries
  void add_to_entries(int cell_index) override
  {
    add_to_entries(cell_index, 1);
  }

  /// increment cell entries
  void add_to_entries(int cell_index, int value) override;

  /// increment left hand side matrix
  void add_to_lhs(int cell_index, int i, int j, float value) override;

  /// increment right hand side column
  void add_to_rhs(int cell_index, int i, float value) override;

  /// add content from other container
  bool add(const TpcSpaceChargeMatrixContainer& other) override;

  //@}

 private:
  /// boundary check
  bool bound_check(int cell_index) const;

  /// boundary check
  bool bound_check(int cell_index, int i) const;

  /// boundary check
  bool bound_check(int cell_index, int i, int j) const;

  /// map matrix index to flat array
  int get_flat_index(int i, int j) const;

  ///@name grid size
  //@{
  int m_pbins = 36; // only for extrapolation
  int m_rbins = 16;
  int m_zbins = 80;
  //@}

  //! number of coordinates
  static constexpr int m_ncoord = 3;

  /// internal matrix representation
  /**
   * Since matrices are symetric, one just needs to store ncoords*(ncoords+1)/2 values
   */
  using matrix_t = std::array<float, m_ncoord * m_ncoord>;
  using column_t = std::array<float, m_ncoord>;

  /// left hand side matrices for distortion inversions
  std::vector<matrix_t> m_lhs;

  /// right hand side matrices for distortion inversions
  std::vector<column_t> m_rhs;

  /// keep track of how many entries are used per cells
  std::vector<int> m_entries;

  ClassDefOverride(TpcSpaceChargeMatrixContainer2D, 1)
};

#endif
