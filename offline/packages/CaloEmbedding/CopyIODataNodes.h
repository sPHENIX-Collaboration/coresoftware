// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef COPYIODATANODES_H
#define COPYIODATANODES_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class CopyIODataNodes : public SubsysReco
{
 public:
  CopyIODataNodes(const std::string &name = "CopyIODataNodes");

  ~CopyIODataNodes() override = default;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  void CopyCentralityInfo(bool flag = true) { m_CopyCentralityInfoFlag = flag; }
  void CopyEventHeader(bool flag = true) { m_CopyEventHeaderFlag = flag; }
  void CopyGlobalVertexMap(bool flag = true) { m_CopyGlobalVertexMapFlag = flag; }
  void CopyMinimumBiasInfo(bool flag = true) { m_CopyMinimumBiasInfoFlag = flag; }
  void CopyRunHeader(bool flag = true) { m_CopyRunHeaderFlag = flag; }
  void CopySyncObject(bool flag = true) { m_CopySyncObjectFlag = flag; }

 private:
  void CreateCentralityInfo(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode);
  void CopyCentralityInfo(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode);

  void CreateEventHeader(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode);
  void CopyEventHeader(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode);

  void CreateGlobalVertexMap(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode);
  void CopyGlobalVertexMap(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode);

  void CreateMinimumBiasInfo(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode);
  void CopyMinimumBiasInfo(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode);

  void CopyRunHeader(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode);

  void CreateSyncObject(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode);
  void CopySyncObject(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode);

  bool m_CopyCentralityInfoFlag = true;
  bool m_CopyEventHeaderFlag = true;
  bool m_CopyGlobalVertexMapFlag = true;
  bool m_CopyMinimumBiasInfoFlag = true;
  bool m_CopyRunHeaderFlag = true;
  bool m_CopySyncObjectFlag = true;
};

#endif  // COPYIODATANODES_H
