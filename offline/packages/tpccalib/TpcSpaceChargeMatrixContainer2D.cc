
/**
 * @file tpccalib/TpcSpaceChargeMatrixContainer2D.cc
 * @author Xudong Yu
 * @date Februray 2025
 * @brief Contains matrices needed for space charge trackbase reconstruction 2D
 */

#include "TpcSpaceChargeMatrixContainer2D.h"

//___________________________________________________________
TpcSpaceChargeMatrixContainer2D::TpcSpaceChargeMatrixContainer2D()
{
  // reset all matrix arrays
  TpcSpaceChargeMatrixContainer2D::Reset();
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainer2D::identify(std::ostream& out) const
{
  out << "TpcSpaceChargeMatrixContainer2D" << std::endl;
  out << "  pbins: " << m_pbins << std::endl;
  out << "  rbins: " << m_rbins << std::endl;
  out << "  zbins: " << m_zbins << std::endl;
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainer2D::get_grid_dimensions(int& pbins, int& rbins, int&  zbins) const
{
  pbins = m_pbins;
  rbins = m_rbins;
  zbins = m_zbins;
}

//___________________________________________________________
int TpcSpaceChargeMatrixContainer2D::get_grid_size() const
{
  return m_rbins * m_zbins;
}

//___________________________________________________________
int TpcSpaceChargeMatrixContainer2D::get_cell_index(int ir, int iz) const
{
  if (ir < 0 || ir >= m_rbins)
  {
    return -1;
  }
  if (iz < 0 || iz >= m_zbins)
  {
    return -1;
  }
  return iz + m_zbins * ir;
}

//___________________________________________________________
int TpcSpaceChargeMatrixContainer2D::get_entries(int cell_index) const
{
  // bound check
  if (!bound_check(cell_index))
  {
    return 0;
  }
  return m_entries[cell_index];
}

//___________________________________________________________
float TpcSpaceChargeMatrixContainer2D::get_lhs(int cell_index, int i, int j) const
{
  // bound check
  if (!bound_check(cell_index, i, j))
  {
    return 0;
  }
  return m_lhs[cell_index][get_flat_index(i, j)];
}

//___________________________________________________________
float TpcSpaceChargeMatrixContainer2D::get_rhs(int cell_index, int i) const
{
  // bound check
  if (!bound_check(cell_index, i))
  {
    return 0;
  }
  return m_rhs[cell_index][i];
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainer2D::Reset()
{
  // reset total number of bins
  const int totalbins = m_rbins * m_zbins;

  // reset arrays
  m_entries = std::vector<int>(totalbins, 0);
  m_lhs = std::vector<matrix_t>(totalbins, {{}});
  m_rhs = std::vector<column_t>(totalbins, {{}});
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainer2D::set_grid_dimensions(int pbins, int rbins, int zbins)
{
  m_pbins = pbins;
  m_rbins = rbins;
  m_zbins = zbins;
  Reset();
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainer2D::add_to_entries(int cell_index, int value)
{
  if (bound_check(cell_index))
  {
    m_entries[cell_index] += value;
  }
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainer2D::add_to_lhs(int cell_index, int i, int j, float value)
{
  if (bound_check(cell_index, i, j))
  {
    m_lhs[cell_index][get_flat_index(i, j)] += value;
  }
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainer2D::add_to_rhs(int cell_index, int i, float value)
{
  if (bound_check(cell_index, i))
  {
    m_rhs[cell_index][i] += value;
  }
}

//___________________________________________________________
bool TpcSpaceChargeMatrixContainer2D::add(const TpcSpaceChargeMatrixContainer& other)
{
  // check dimensions
  int pbins = 0;
  int rbins = 0;
  int zbins = 0;
  other.get_grid_dimensions(pbins, rbins, zbins);
  if ((m_pbins != pbins) || (m_rbins != rbins) || (m_zbins != zbins))
  {
    std::cout << "TpcSpaceChargeMatrixContainer2D::add - inconsistent grid sizes" << std::endl;
    return false;
  }

  // increment cell entries
  for (size_t cell_index = 0; cell_index < m_lhs.size(); ++cell_index)
  {
    add_to_entries(cell_index, other.get_entries(cell_index));
  }

  // increment left hand side matrices
  for (size_t cell_index = 0; cell_index < m_lhs.size(); ++cell_index)
  {
    for (int i = 0; i < m_ncoord; ++i)
    {
      for (int j = 0; j < m_ncoord; ++j)
      {
        add_to_lhs(cell_index, i, j, other.get_lhs(cell_index, i, j));
      }
    }
  }

  // increment right hand side matrices
  for (size_t cell_index = 0; cell_index < m_lhs.size(); ++cell_index)
  {
    for (int i = 0; i < m_ncoord; ++i)
    {
      add_to_rhs(cell_index, i, other.get_rhs(cell_index, i));
    }
  }

  return true;
}

//___________________________________________________________
bool TpcSpaceChargeMatrixContainer2D::bound_check(int cell_index) const
{
  if (cell_index < 0 || cell_index >= (int) m_rhs.size())
  {
    return false;
  }
  return true;
}

//___________________________________________________________
bool TpcSpaceChargeMatrixContainer2D::bound_check(int cell_index, int i) const
{
  if (cell_index < 0 || cell_index >= (int) m_rhs.size())
  {
    return false;
  }
  if (i < 0 || i >= m_ncoord)
  {
    return false;
  }
  return true;
}

//___________________________________________________________
bool TpcSpaceChargeMatrixContainer2D::bound_check(int cell_index, int i, int j) const
{
  if (cell_index < 0 || cell_index >= (int) m_lhs.size())
  {
    return false;
  }
  if (i < 0 || i >= m_ncoord)
  {
    return false;
  }
  if (j < 0 || j >= m_ncoord)
  {
    return false;
  }
  return true;
}

//___________________________________________________________
int TpcSpaceChargeMatrixContainer2D::get_flat_index(int i, int j) const
{
  return j + i * m_ncoord;
}
