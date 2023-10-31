#include "DumpMbdOut.h"

#include <mbd/MbdOut.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<MbdOut>;

DumpMbdOut::DumpMbdOut(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpMbdOut::process_Node(PHNode *myNode)
{
  MbdOut *bbcout = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    bbcout = thisNode->getData();
  }
  if (bbcout && bbcout->isValid())
  {
    *fout << "MbdOut->get_bz: " << bbcout->get_zvtx() << std::endl;
    *fout << "MbdOut->get_bzerr: " << bbcout->get_zvtxerr() << std::endl;
    *fout << "MbdOut->get_t0: " << bbcout->get_t0() << std::endl;
    *fout << "MbdOut->get_t0err: " << bbcout->get_t0err() << std::endl;
    *fout << "MbdOut->get_bns: " << bbcout->get_npmt(0) << std::endl;
    *fout << "MbdOut->get_bnn: " << bbcout->get_npmt(1) << std::endl;
    *fout << "MbdOut->get_bqs: " << bbcout->get_q(0) << std::endl;
    *fout << "MbdOut->get_bqn: " << bbcout->get_q(1) << std::endl;
    *fout << "MbdOut->get_bts: " << bbcout->get_time(0) << std::endl;
    *fout << "MbdOut->get_btn: " << bbcout->get_time(1) << std::endl;
  }
  return 0;
}
