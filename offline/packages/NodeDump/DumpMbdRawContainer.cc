#include "DumpMbdRawContainer.h"

#include <mbd/MbdRawContainer.h>
#include <mbd/MbdRawHit.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<MbdRawContainer>;

DumpMbdRawContainer::DumpMbdRawContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpMbdRawContainer::process_Node(PHNode *myNode)
{
  MbdRawContainer *mbdraw = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);  // NOLINT(cppcoreguidelines-pro-type-static-cast-downcast)
  if (thisNode)
  {
    mbdraw = thisNode->getData();
  }
  if (mbdraw && mbdraw->isValid())
  {
    int npmt = mbdraw->get_npmt();
    *fout << "MbdRawContainer->get_evt(): " << mbdraw->get_evt()<< std::endl;
    *fout << "MbdRawContainer->get_clock(): " << mbdraw->get_clock() << std::endl;
    *fout << "MbdRawContainer->get_femclock(): " << mbdraw->get_femclock() << std::endl;
    *fout << "MbdRawContainer->get_npmt(): " << npmt << std::endl;
    *fout << "MbdRawContainer->get_badtdc(): " << mbdraw->get_badtdc() << std::endl;
    
    for (int j = 0; j < npmt; j++)
    {
      MbdRawHit *rawhit = mbdraw->get_pmt(j);
      *fout << "MbdRawHit[" << j << "]->get_pmt(): " << rawhit->get_pmt() << std::endl;
      *fout << "MbdRawHit[" << j << "]->get_adc(): " << rawhit->get_adc() << std::endl;
      *fout << "MbdRawHit[" << j << "]->get_ttdc(): " << rawhit->get_ttdc() << std::endl;
      *fout << "MbdRawHit[" << j << "]->get_qtdc(): " << rawhit->get_qtdc() << std::endl;
      *fout << "MbdRawHit[" << j << "]->get_chi2ndf(): " << rawhit->get_chi2ndf() << std::endl;
      *fout << "MbdRawHit[" << j << "]->get_fitinfo(): " << rawhit->get_fitinfo() << std::endl;
      *fout << "MbdRawHit[" << j << "]->get_badtdc(): " << rawhit->get_badtdc() << std::endl;
      *fout << "MbdRawHit[" << j << "]->get_fitstat(): " << rawhit->get_fitstat() << std::endl;
    }
  }
  return 0;
}
