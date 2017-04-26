// TPC CHANNEL class
// Stores TPC detection pad object
// Author: Carlos Perez
#ifndef __TPCCHANNEL_H__
#define __TPCCHANNEL_H__

class TPCChannel : public vChannel {
 public:
  TPCChannel() {}
  virtual ~TPCChannel() {}

  UChar_t GetModuleId() {return fID1;}
  void SetModuleId(UChar_t id) {fID1=id;}
  void SetModuleId(int ra, int ph, int zz) {fID1=36*zz + 12*ra + ph;}
  UChar_t GetModR() {return fID1/12;}
  UChar_t GetModP() {return fID1%12;}
  UChar_t GetModZ() {return fID1/36;}

  UShort_t GetPRCId() {return fID2;}
  void SetPRCId(UShort_t id) {fID2=id;}
  void SetPRC(int row, int col) {fID2=8*col + row;}
  UShort_t GetPRow() {return fID2%8;}
  UShort_t GetPCol() {return fID2/8;}

 protected:
};

#endif
