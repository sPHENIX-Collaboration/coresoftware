#ifndef FUN4ALLPACKET_LL1PACKETCONTAINERV1_H
#define FUN4ALLPACKET_LL1PACKETCONTAINERV1_H

#include "LL1PacketContainer.h"

class LL1Packet;
class TClonesArray;

class LL1PacketContainerv1 : public LL1PacketContainer
{
 public:
  LL1PacketContainerv1();
  ~LL1PacketContainerv1() override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  LL1Packet *AddPacket() override;
  LL1Packet *AddPacket(LL1Packet *ll1packet) override;
  unsigned int get_npackets() override;
  LL1Packet *getPacket(unsigned int index) override;
  LL1Packet *getPacketbyId(int id) override;
  void setEvtSequence(const int i) override {eventno = i;}
  int getEvtSequence() const override {return eventno;}
  void setStatus(const unsigned int ui) override {status = ui;}
  unsigned int getStatus() const override {return status;}

 private:
  TClonesArray *LL1PacketsTCArray{nullptr};
  int eventno{std::numeric_limits<int>::min()};
  unsigned int status {0};

  ClassDefOverride(LL1PacketContainerv1, 2)
};

#endif
