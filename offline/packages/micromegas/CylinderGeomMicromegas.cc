/*!
 * \file CylinderGeomMicromegas.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "CylinderGeomMicromegas.h"

void CylinderGeomMicromegas::identify( std::ostream& out ) const
{
  out << "CylinderGeomMicromegas" << std::endl;
  out << "layer: " << m_layer << std::endl;
  out << "segmentation_type: " << (m_segmentation_type == MicromegasDefs::SegmentationType::SEGMENTATION_PHI ? "SEGMENTATION_PHI":"SEGMENTATION_Z") << std::endl;
  out << "radius: " << m_radius << "cm" << std::endl;
  out << "thickness: " << m_thickness << "cm" << std::endl;
  out << "zmin: " << m_zmin << "cm" << std::endl;
  out << "zmax: " << m_zmax << "cm" << std::endl;
  out << std::endl;
}
