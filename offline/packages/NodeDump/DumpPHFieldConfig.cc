#include "DumpPHFieldConfig.h"

#include "DumpObject.h"

#include <phfield/PHFieldConfig.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<PHFieldConfig>;

DumpPHFieldConfig::DumpPHFieldConfig(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHFieldConfig::process_Node(PHNode *myNode)
{
  PHFieldConfig *phfieldconfig = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    phfieldconfig = thisNode->getData();
  }
  if (phfieldconfig)
  {
    *fout << "PHFieldConfig->isValid(): " << phfieldconfig->isValid() << std::endl;
    if (phfieldconfig->isValid())
    {
      *fout << "get_field_config(): " << phfieldconfig->get_field_config() << std::endl;
      *fout << "get_filename(): " << phfieldconfig->get_filename() << std::endl;
      *fout << "get_magfield_rescale(): " << phfieldconfig->get_magfield_rescale() << std::endl;
      *fout << "get_field_mag_x(): " << phfieldconfig->get_field_mag_x() << std::endl;
      *fout << "get_field_mag_y(): " << phfieldconfig->get_field_mag_y() << std::endl;
      *fout << "get_field_mag_z(): " << phfieldconfig->get_field_mag_z() << std::endl;
    }
  }
  return 0;
}
