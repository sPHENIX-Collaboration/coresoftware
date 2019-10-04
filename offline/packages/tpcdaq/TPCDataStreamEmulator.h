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

#include <vector>

class PHCompositeNode;
class Fun4AllHistoManager;
class TH1;
class TH2;

/*!
 * \brief TPCDataStreamEmulator
 */
class TPCDataStreamEmulator : public SubsysReco
{
 public:
  TPCDataStreamEmulator(
      unsigned int minLayer,
      unsigned int m_maxLayer,
      const std::string &outputfilenamebase = "TPCDataStreamEmulator");
  virtual ~TPCDataStreamEmulator();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void maxLayer(int maxLayer)
  {
    m_maxLayer = maxLayer;
  }

  void minLayer(int minLayer)
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
#if !defined(__CINT__) || defined(__CLING__)

  Fun4AllHistoManager *getHistoManager();
  int writeWavelet(int layer, int side, int phibin, int hittime, const std::vector<unsigned int> &wavelet);

  bool m_saveDataStreamFile;

  std::string m_outputFileNameBase;

  int m_minLayer;
  int m_maxLayer;

  int m_evtCounter;

  double m_vertexZAcceptanceCut;
  double m_etaAcceptanceCut;

  // histograms
  TH1 *m_hDataSize;
  TH1 *m_hWavelet;
  TH1 *m_hNChEta;
  TH2 *m_hLayerWaveletSize;
  TH2 *m_hLayerHit;
  TH2 *m_hLayerZBinHit;
  TH2 *m_hLayerZBinADC;
  TH2 *m_hLayerDataSize;
  TH2 *m_hLayerSumHit;
  TH2 *m_hLayerSumDataSize;

#endif  // #if !defined(__CINT__) || defined(__CLING__)
};

#endif /* TPCDATASTREAMEMULATOR_H_ */
