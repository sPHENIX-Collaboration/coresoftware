#include "DumpPacket.h"

#include <ffarawobjects/OfflinePacket.h>

DumpPacket::~DumpPacket()
{
  for (auto &iter: m_PacketDumpFile)
  {
    iter.second->close();
  }
  return;
}

void DumpPacket::ddumppacket(OfflinePacket *pkt)
{
  int packetid = pkt->getIdentifier();
  if (m_PacketDumpFile.find(packetid) == m_PacketDumpFile.end())
  {
    std::string fname = "offlinepacket_" + std::to_string(packetid) + ".ddump";
    std::ofstream *dumpfile = new std::ofstream(fname);
    //dumpfile.open(fname);
      m_PacketDumpFile.insert(std::make_pair(packetid,dumpfile));
      m_PacketDumpCounter.insert(std::make_pair(packetid,m_ddump_flag));
  }
  if (m_PacketDumpCounter[packetid] != 0)
  {
  pkt->dump(*m_PacketDumpFile[packetid]);
    m_PacketDumpCounter[packetid]--;
  }
  return;
}
