
/**
 * @file tpccalib/TpcSpaceChargeMatrixContainerv2.cc
 * @author Hugo Pereira Da Costa
 * @date June 2018
 * @brief Contains matrices needed for space charge trackbase reconstruction
 */

#include "TpcSpaceChargeMatrixContainerv2.h"

//___________________________________________________________
TpcSpaceChargeMatrixContainerv2::TpcSpaceChargeMatrixContainerv2()
{
  // reset all matrix arrays
  TpcSpaceChargeMatrixContainerv2::Reset();
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::identify(std::ostream& out) const
{
  out << "TpcSpaceChargeMatrixContainerv2" << std::endl;
  out << "  phibins: " << m_phibins << std::endl;
  out << "  rbins: " << m_rbins << std::endl;
  out << "  zbins: " << m_zbins << std::endl;
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::get_grid_dimensions(int& phibins, int& rbins, int& zbins) const
{
  phibins = m_phibins;
  rbins = m_rbins;
  zbins = m_zbins;
}

//___________________________________________________________
int TpcSpaceChargeMatrixContainerv2::get_grid_size() const
{
  return m_phibins * m_rbins * m_zbins;
}

//___________________________________________________________
int TpcSpaceChargeMatrixContainerv2::get_cell_index(int iphi, int ir, int iz) const
{
  if (iphi < 0 || iphi >= m_phibins)
  {
    return -1;
  }
  if (ir < 0 || ir >= m_rbins)
  {
    return -1;
  }
  if (iz < 0 || iz >= m_zbins)
  {
    return -1;
  }
  return iz + m_zbins * (ir + m_rbins * iphi);
}

//___________________________________________________________
int TpcSpaceChargeMatrixContainerv2::get_entries(int cell_index) const
{
  // bound check
  if (!bound_check(cell_index))
  {
    return 0;
  }
  return m_entries[cell_index];
}

//___________________________________________________________
float TpcSpaceChargeMatrixContainerv2::get_lhs(int cell_index, int i, int j) const
{
  // bound check
  if (!bound_check(cell_index, i, j))
  {
    return 0;
  }
  return m_lhs[cell_index][get_flat_index(i, j)];
}

//___________________________________________________________
float TpcSpaceChargeMatrixContainerv2::get_rhs(int cell_index, int i) const
{
  // bound check
  if (!bound_check(cell_index, i))
  {
    return 0;
  }
  return m_rhs[cell_index][i];
}

//___________________________________________________________
float TpcSpaceChargeMatrixContainerv2::get_lhs_rphi(int cell_index, int i, int j) const
{
  // bound check
  if (!bound_check_reduced(cell_index, i, j))
  {
    return 0;
  }
  return m_lhs_rphi[cell_index][get_flat_index_reduced(i, j)];
}

//___________________________________________________________
float TpcSpaceChargeMatrixContainerv2::get_rhs_rphi(int cell_index, int i) const
{
  // bound check
  if (!bound_check_reduced(cell_index, i))
  {
    return 0;
  }
  return m_rhs_rphi[cell_index][i];
}

//___________________________________________________________
float TpcSpaceChargeMatrixContainerv2::get_lhs_z(int cell_index, int i, int j) const
{
  // bound check
  if (!bound_check_reduced(cell_index, i, j))
  {
    return 0;
  }
  return m_lhs_z[cell_index][get_flat_index_reduced(i, j)];
}

//___________________________________________________________
float TpcSpaceChargeMatrixContainerv2::get_rhs_z(int cell_index, int i) const
{
  // bound check
  if (!bound_check_reduced(cell_index, i))
  {
    return 0;
  }
  return m_rhs_z[cell_index][i];
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::Reset()
{
  // reset total number of bins
  const int totalbins = m_phibins * m_rbins * m_zbins;

  // reset arrays
  m_entries = std::vector<int>(totalbins, 0);
  m_lhs = std::vector<matrix_t>(totalbins, {{}});
  m_rhs = std::vector<column_t>(totalbins, {{}});

  // rphi reduced arrays
  m_lhs_rphi = std::vector<reduced_matrix_t>(totalbins, {{}});
  m_rhs_rphi = std::vector<reduced_column_t>(totalbins, {{}});

  // z reduced arrays
  m_lhs_z = std::vector<reduced_matrix_t>(totalbins, {{}});
  m_rhs_z = std::vector<reduced_column_t>(totalbins, {{}});
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::set_grid_dimensions(int phibins, int rbins, int zbins)
{
  m_phibins = phibins;
  m_rbins = rbins;
  m_zbins = zbins;
  Reset();
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::add_to_entries(int cell_index, int value)
{
  if (bound_check(cell_index))
  {
    m_entries[cell_index] += value;
  }
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::add_to_lhs(int cell_index, int i, int j, float value)
{
  if (bound_check(cell_index, i, j))
  {
    m_lhs[cell_index][get_flat_index(i, j)] += value;
  }
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::add_to_rhs(int cell_index, int i, float value)
{
  if (bound_check(cell_index, i))
  {
    m_rhs[cell_index][i] += value;
  }
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::add_to_lhs_rphi(int cell_index, int i, int j, float value)
{
  if (bound_check_reduced(cell_index, i, j))
  {
    m_lhs_rphi[cell_index][get_flat_index_reduced(i, j)] += value;
  }
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::add_to_rhs_rphi(int cell_index, int i, float value)
{
  if (bound_check_reduced(cell_index, i))
  {
    m_rhs_rphi[cell_index][i] += value;
  }
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::add_to_lhs_z(int cell_index, int i, int j, float value)
{
  if (bound_check_reduced(cell_index, i, j))
  {
    m_lhs_z[cell_index][get_flat_index_reduced(i, j)] += value;
  }
}

//___________________________________________________________
void TpcSpaceChargeMatrixContainerv2::add_to_rhs_z(int cell_index, int i, float value)
{
  if (bound_check_reduced(cell_index, i))
  {
    m_rhs_z[cell_index][i] += value;
  }
}

//___________________________________________________________
bool TpcSpaceChargeMatrixContainerv2::add(const TpcSpaceChargeMatrixContainer& other)
{
  // check dimensions
  int phibins = 0;
  int rbins = 0;
  int zbins = 0;
  other.get_grid_dimensions(phibins, rbins, zbins);
  if ((m_phibins != phibins) || (m_rbins != rbins) || (m_zbins != zbins))
  {
    std::cout << "TpcSpaceChargeMatrixContainerv2::add - inconsistent grid sizes" << std::endl;
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
  for (size_t cell_index = 0; cell_index < m_rhs.size(); ++cell_index)
  {
    for (int i = 0; i < m_ncoord; ++i)
    {
      add_to_rhs(cell_index, i, other.get_rhs(cell_index, i));
    }
  }

  // increment reduced rphi left hand side matrices
  for (size_t cell_index = 0; cell_index < m_lhs_rphi.size(); ++cell_index)
  {
    for (int i = 0; i < m_ncoord_reduced; ++i)
    {
      for (int j = 0; j < m_ncoord_reduced; ++j)
      {
        add_to_lhs_rphi(cell_index, i, j, other.get_lhs_rphi(cell_index, i, j));
      }
    }
  }

  // increment reduced rphi right hand side matrices
  for (size_t cell_index = 0; cell_index < m_rhs_rphi.size(); ++cell_index)
  {
    for (int i = 0; i < m_ncoord_reduced; ++i)
    {
      add_to_rhs_rphi(cell_index, i, other.get_rhs_rphi(cell_index, i));
    }
  }

  // increment reduced z left hand side matrices
  for (size_t cell_index = 0; cell_index < m_lhs_z.size(); ++cell_index)
  {
    for (int i = 0; i < m_ncoord_reduced; ++i)
    {
      for (int j = 0; j < m_ncoord_reduced; ++j)
      {
        add_to_lhs_z(cell_index, i, j, other.get_lhs_z(cell_index, i, j));
      }
    }
  }

  // increment reduced z right hand side matrices
  for (size_t cell_index = 0; cell_index < m_rhs_z.size(); ++cell_index)
  {
    for (int i = 0; i < m_ncoord_reduced; ++i)
    {
      add_to_rhs_z(cell_index, i, other.get_rhs_z(cell_index, i));
    }
  }

  return true;
}

//___________________________________________________________
bool TpcSpaceChargeMatrixContainerv2::bound_check(int cell_index) const
{
  if (cell_index < 0 || cell_index >= (int) m_rhs.size())
  {
    return false;
  }
  return true;
}

//___________________________________________________________
bool TpcSpaceChargeMatrixContainerv2::bound_check(int cell_index, int i) const
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
bool TpcSpaceChargeMatrixContainerv2::bound_check(int cell_index, int i, int j) const
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
bool TpcSpaceChargeMatrixContainerv2::bound_check_reduced(int cell_index, int i) const
{
  if (cell_index < 0 || cell_index >= (int) m_rhs_rphi.size())
  {
    return false;
  }
  if (i < 0 || i >= m_ncoord_reduced)
  {
    return false;
  }
  return true;
}

//___________________________________________________________
bool TpcSpaceChargeMatrixContainerv2::bound_check_reduced(int cell_index, int i, int j) const
{
  if (cell_index < 0 || cell_index >= (int) m_lhs_rphi.size())
  {
    return false;
  }
  if (i < 0 || i >= m_ncoord_reduced)
  {
    return false;
  }
  if (j < 0 || j >= m_ncoord_reduced)
  {
    return false;
  }
  return true;
}

//___________________________________________________________
int TpcSpaceChargeMatrixContainerv2::get_flat_index(int i, int j) const
{
  return j + i * m_ncoord;
}

//___________________________________________________________
int TpcSpaceChargeMatrixContainerv2::get_flat_index_reduced(int i, int j) const
{
  return j + i * m_ncoord_reduced;
}
