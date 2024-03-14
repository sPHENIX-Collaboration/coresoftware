#include <fstream>
#include <map>

class OfflinePacket;

#ifndef FFARAWMODULES_DUMPPACKET_H
#define FFARAWMODULES_DUMPPACKET_H

class DumpPacket
{
public:
  DumpPacket() = default;
  virtual ~DumpPacket() = default;
  virtual void ddumppacket(OfflinePacket *pkt);
  virtual void enable_ddump(int b = 1) {m_ddump_flag = b;}
  virtual bool ddump_enabled() const {return m_ddump_flag;}

private:
  int m_ddump_flag {0};
  std::map<int, std::ofstream *> m_PacketDumpFile;
  std::map<int, int> m_PacketDumpCounter;
};

#endif
