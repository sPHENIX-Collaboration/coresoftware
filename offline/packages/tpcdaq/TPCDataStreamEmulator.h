// $Id: $

/*!
 * \file TPCDataStreamEmulator.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef TPCDATASTREAMEMULATOR_H_
#define TPCDATASTREAMEMULATOR_H_

#include <fun4all/SubsysReco.h>

class PHCompositeNode;
class Fun4AllHistoManager;

/*!
 * \brief TPCDataStreamEmulator
 */
class TPCDataStreamEmulator : public SubsysReco
{
 public:
  TPCDataStreamEmulator(
      unsigned int minLayer,
      unsigned int m_maxLayer,
      const std::string &outputfilenamebase = "TPCIntegratedCharge.root");
  virtual ~TPCDataStreamEmulator();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void maxLayer(unsigned int maxLayer)
  {
    m_maxLayer = maxLayer;
  }

  void minLayer(unsigned int minLayer)
  {
    m_minLayer = minLayer;
  }

  void outputFileNameBase(const std::string &outputFileNameBase)
  {
    m_outputFileNameBase = outputFileNameBase;
  }

  void saveDataStreamFile(bool saveDataStreamFile)
  {
    m_saveDataStreamFile = saveDataStreamFile;
  }

 private:
#ifndef __CINT__

  Fun4AllHistoManager *getHistoManager();
  void writeWavelet(int layer, int side, int phibin, int hittime, const std::vector<unsigned int>  &last_wavelet);

  bool m_saveDataStreamFile;

  std::string m_outputFileNameBase;

  unsigned int m_minLayer;
  unsigned int m_maxLayer;

  int m_evtCounter;

#endif  // #ifndef __CINT__
};

#endif /* TPCDATASTREAMEMULATOR_H_ */
