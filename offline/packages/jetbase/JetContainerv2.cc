#include "JetContainerv2.h"
#include "Jetv3.h"

JetContainerv2::JetContainerv2()
  : JetContainerv1("Jetv3")
{
}

void JetContainerv2::identify(std::ostream& os) const
{
  os << "JetContainerv2 (Jetv3 jets with corrected pT): size = " << size() << std::endl;
  if (!get_vertex_type().empty())
  {
    os << "  reconstructed with vertex type " << get_vertex_type() << " at z = " << get_vertex_z() << std::endl;
  }
  os << "  Contains jets with the following properties:" << std::endl;
  print_property_types(os);
}
