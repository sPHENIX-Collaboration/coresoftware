#ifndef INTT_COMBINEDRAWDATADECODER_H
#define INTT_COMBINEDRAWDATADECODER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class InttCombinedRawDataDecoder : public SubsysReco
{
 public:
  InttCombinedRawDataDecoder(std::string const& name = "InttCombinedRawDataDecoder");

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

 private:
  std::string m_InttRawNodeName = "INTTRAWHIT";
};

#endif  // INTT_COMBINEDRAWDATADECODER_H
