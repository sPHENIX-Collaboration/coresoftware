#ifndef FUN4ALLRAW_SINGLEGL1TRIGGEREDINPUT_H
#define FUN4ALLRAW_SINGLEGL1TRIGGEREDINPUT_H

#include "SingleTriggeredInput.h"

#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <vector>

class OfflinePacket;
class Packet;
class PHCompositeNode;

class SingleGl1TriggeredInput : public SingleTriggeredInput
{
 public:
  explicit SingleGl1TriggeredInput(const std::string &name);
  ~SingleGl1TriggeredInput() override = default;
  void FillPool() override;
  // void CleanupUsedPackets(const int eventno);
  // void ClearCurrentEvent();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNodes(Event *evt) override;
  uint64_t GetClock(Event *evt, int pid) override;
  int GetCurrentPacketNumber() const { return m_PacketNumber; }
  int GetLastPacketNumber() const { return m_LastPacketNumber; }
  const std::array<int, pooldepth>& GetGl1SkipArray() const { return m_Gl1SkipPerIndex; }
  const std::array<uint64_t, pooldepth>& GetPacketNumbers() const { return m_Gl1PacketNumbers; }
  int ReadEvent() override;
  void SetPacketNumbers(int last, int current)
  {
    m_LastPacketNumber = last;
    m_PacketNumber = current;
  }
  void SetGl1SkipAtIndex(size_t index, int value)
  {
    if (index < pooldepth)
    {
      m_Gl1SkipPerIndex[index] = value;
    }
  }
  void SetGl1PacketNumber(size_t index, uint64_t value)
  {
    if (index < pooldepth)
    {
      m_Gl1PacketNumbers[index] = value;
    }
  }


 protected:
  int m_LastPacketNumber{0};
  int m_PacketNumber{0};
  std::array<int, pooldepth> m_Gl1SkipPerIndex{};
  std::array<uint64_t, pooldepth> m_Gl1PacketNumbers{};

 private:
};

#endif
