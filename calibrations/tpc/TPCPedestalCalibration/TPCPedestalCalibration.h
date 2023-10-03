// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCPedestalCalibration_H
#define TPCPedestalCalibration_H

#include <fun4all/SubsysReco.h>
#include <cdbobjects/CDBTTree.h>
#include <sphenixnpc/CDBUtils.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;

class TPCPedestalCalibration : public SubsysReco
{
 public:
  explicit TPCPedestalCalibration(const std::string &name = "TPCPedestalCalibration.root");

  ~TPCPedestalCalibration() override {}

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int EndRun(const int runnumber) override;

  void CDBInsert();

  int End(PHCompositeNode *topNode) override; 

  void AddPacket(int packet)
  {
    m_packets.push_back(packet);
  }

  void SetSector(int sectorNum)
  {
    m_sector = sectorNum;
  }

  void InsertCDBTTree(std::string username) // username for CDB record of who uploaded the file
  {
    m_writeToCDB = true;
    m_username = username; 
  }

  void NewCDBGlobalTag(std::string username)
  {
    CDBUtils *uti = new CDBUtils();
    uti->createGlobalTag(username); 
  }
 
 protected:
  //! which packet to decode
  std::vector<int> m_packets{1001};

 private:
  std::string m_fname;
  bool m_writeToCDB;
  CDBTTree * m_cdbttree = nullptr;

  int m_BCO = 0;
  int m_packet = 0;
  int m_nWaveormInFrame = 0;
  int m_nSamples = 0;
  int m_fee = 0;
  int m_Channel = 0;
  std::vector<unsigned short> m_adcSamples; 

  float m_aveADCFeeChannel[26][256];
  float m_stdADCFeeChannel[26][256];
  float m_countsADCFeeChannel[26][256];
  int m_aliveArrayFeeChannel[26][256];

  std::string m_username = "test";
  bool m_firstBCO = true;
  int m_isAlive = 1;
  float m_pedMean = 0;
  float m_pedStd = 0;
  int m_sector = 0;
  int m_outFEE = 99;
  int m_chan = 999;
  int m_module = 9;
  int m_slot = 99;
 
  int mod_arr[26] = {2,2,1,1,1,3,3,3,3,3,3,2,2,1,2,2,1,1,2,2,3,3,3,3,3,3};
  int slot_arr[26] = {5,6,1,3,2,12,10,11,9,8,7,1,2,4,8,7,6,5,4,3,1,3,2,4,6,5};
};

#endif // TPCPedestalCalibration_H
