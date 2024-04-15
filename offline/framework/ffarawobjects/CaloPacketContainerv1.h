#ifndef FUN4ALLPACKET_CALOPACKETCONTAINERV1_H
#define FUN4ALLPACKET_CALOPACKETCONTAINERV1_H

#include "CaloPacketContainer.h"

class CaloPacket;
class TClonesArray;

class CaloPacketContainerv1 : public CaloPacketContainer
{
 public:
  CaloPacketContainerv1();
  ~CaloPacketContainerv1() override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  CaloPacket *AddPacket() override;
  CaloPacket *AddPacket(CaloPacket *calopacket) override;
  unsigned int get_npackets() override;
  CaloPacket *getPacket(unsigned int index) override;
  CaloPacket *getPacketbyId(int id) override;

 private:
  TClonesArray *CaloPacketsTCArray{nullptr};

  ClassDefOverride(CaloPacketContainerv1, 1)
};

#endif
