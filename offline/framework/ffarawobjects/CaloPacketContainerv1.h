#ifndef FUN4ALLPACKET_CALOPACKETCONTAINERV1_H
#define FUN4ALLPACKET_CALOPACKETCONTAINERV1_H

#include "CaloPacketContainer.h"

#include <limits>

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
  void setEvtSequence(const int i) override {eventno = i;}
  int getEvtSequence() const override {return eventno;}
  void setStatus(const unsigned int ui) override {status = ui;}
  unsigned int getStatus() const override {return status;}
  void deletePacketAt(int index) override;
  void deletePacket(CaloPacket *packet) override;

 private:
  TClonesArray *CaloPacketsTCArray{nullptr};
  int eventno{std::numeric_limits<int>::min()};
  unsigned int status {0};

  ClassDefOverride(CaloPacketContainerv1, 2)
};

#endif
