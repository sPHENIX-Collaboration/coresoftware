#ifndef FUN4ALLPACKET_MBDPACKETCONTAINERV1_H
#define FUN4ALLPACKET_MBDPACKETCONTAINERV1_H

#include "MbdPacketContainer.h"

class MbdPacket;
class TClonesArray;

class  MbdPacketContainerv1: public MbdPacketContainer
{
public:
  MbdPacketContainerv1();
  ~MbdPacketContainerv1() override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  MbdPacket *AddPacket() override;
  MbdPacket *AddPacket(MbdPacket *mbdpacket) override;
  unsigned int get_npackets() override;
  MbdPacket *getPacket(unsigned int index) override;

private:
   TClonesArray *MbdPacketsTCArray {nullptr};

  ClassDefOverride(MbdPacketContainerv1,1)
};

#endif
