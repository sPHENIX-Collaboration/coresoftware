#ifndef __PHG4INTDIGITIZER__
#define __PHG4INTTDIGITIZER__

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

#include <assert.h>
#include <float.h>
#include <vector>

class SvtxHitMap;

class PHG4INTTDigitizer : public SubsysReco
{
 public:
  PHG4INTTDigitizer(const std::string &name = "PHG4INTTDigitizer");
  virtual ~PHG4INTTDigitizer() {}

  //! module initialization
  int Init(PHCompositeNode *topNode) { return 0; }

  //! run initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  //! end of process
  int End(PHCompositeNode *topNode);

  void set_adc_scale(const int &layer, const std::vector<double> &userrange)
  {
    if (userrange.size() != nadcbins)
      assert(!"Error: vector in set_fphx_adc_scale(vector) must have eight elements.");

    //sort(userrange.begin(), userrange.end()); // TODO, causes GLIBC error

    std::vector<std::pair<double, double> > vadcrange;
    for (unsigned int irange = 0; irange < userrange.size(); ++irange)
      if (irange == userrange.size() - 1)
        vadcrange.push_back(std::make_pair(userrange[irange], FLT_MAX));
      else
        vadcrange.push_back(std::make_pair(userrange[irange], userrange[irange + 1]));

    _max_fphx_adc.insert(std::make_pair(layer, vadcrange));
  }

 private:
  void CalculateLadderCellADCScale(PHCompositeNode *topNode);

  void DigitizeLadderCells(PHCompositeNode *topNode);
  void PrintHits(PHCompositeNode *topNode);

  // settings
  std::map<int, unsigned int> _max_adc;
  std::map<int, float> _energy_scale;

  // storage
  SvtxHitMap *_hitmap;

  PHTimeServer::timer _timer;  ///< Timer

  const unsigned int nadcbins = 8;
  std::map<int, std::vector<std::pair<double, double> > > _max_fphx_adc;

  unsigned int m_nCells;
  unsigned int m_nDeadCells;
};

#endif
