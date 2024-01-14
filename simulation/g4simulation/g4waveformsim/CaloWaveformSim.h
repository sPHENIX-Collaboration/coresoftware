// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOWAVEFORMSIM_H
#define CALOWAVEFORMSIM_H

#include <fun4all/SubsysReco.h>
#include <caloreco/CaloTowerDefs.h>
#include <calobase/TowerInfoDefs.h>

#include <g4detectors/LightCollectionModel.h> 

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TProfile;
class PHG4Hit;
class PHG4CylinderCellGeom_Spacalv1;
class PHG4CylinderGeom_Spacalv3;
class TRandom3;
class TTree;
class CDBTTree;
class TowerInfoContainer;

class CaloWaveformSim : public SubsysReco
{
public:
    CaloWaveformSim(const std::string &name = "CaloWaveformSim");

    ~CaloWaveformSim() override;

    enum NoiseType{
        NOISE_NONE = 0,
        NOISE_GAUSSIAN = 1,
        NOISE_TREE = 2
    };
    

    int Init(PHCompositeNode *topNode) override;

    /** Called for each event.
        This is where you do the real work.
     */
    int process_event(PHCompositeNode *topNode) override;

    /// Called at the end of all processing.
    int End(PHCompositeNode *topNode) override;

    void set_detector_type(CaloTowerDefs::DetectorSystem dettype)
    {
        m_dettype = dettype;
        return;
    }
    void set_detector(const std::string &detector)
    {
        m_detector = detector;
        return;
    }
    void set_fieldname(const std::string &fieldname)
    {
        m_fieldname = fieldname;
        return;
    }
    void set_calibName(const std::string &calibName)
    {
        m_calibName = calibName;
        m_overrideCalibName = true;
        return;
    }
    void set_overrideCalibName(bool overrideCalibName)
    {
        m_overrideCalibName = overrideCalibName;
        return;
    }
    void set_overrideFieldName(bool overrideFieldName)
    {
        m_overrideFieldName = overrideFieldName;
        return;
    }
    void set_templatefile(const std::string &templatefile)
    {
        m_templatefile = templatefile;
        return;
    }
    void set_nsamples(int _nsamples)
    {
        m_nsamples = _nsamples;
        return;
    }
    void set_noise_type(NoiseType noiseType)
    {
        m_noiseType = noiseType;
        return;
    }
    void set_fixpedestal(int _fixpedestal)
    {
        m_fixpedestal = _fixpedestal;
        return;
    }
    void set_gaussian_noise(int _gaussian_noise)
    {
        m_gaussian_noise = _gaussian_noise;
        return;
    }
    void set_deltaT(float _deltaT)
    {
        m_deltaT = _deltaT;
        return;
    }
    void set_timewidth(float _timewidth)
    {
      m_timeshiftwidth = _timewidth;
      return;
    }
    void set_peakpos(float _peakpos)
    {
      m_peakpos = _peakpos;
      return;
    }
    //for CEMC light yield correction
    LightCollectionModel &get_light_collection_model() { return light_collection_model; }


private:
    CaloTowerDefs::DetectorSystem m_dettype{CaloTowerDefs::CEMC};
    std::string m_detector{"CEMC"};
    std::string m_fieldname{"Femc_datadriven_qm1_correction"};
    std::string m_calibName{"cemc_pi0_twrSlope_v1"};
    std::string m_noisetree_name{"/sphenix/user/shuhangli/noisetree/macro/noiseout_emcalout.root"};
    bool m_overrideCalibName{false};
    bool m_overrideFieldName{false};
    CDBTTree *cdbttree{nullptr};
    std::string m_templatefile{"/sphenix/user/shuhangli/cosmicreco/macro/waveformtemptempohcalcosmic.root"};
    TProfile *h_template{nullptr};
    TTree *m_noisetree{nullptr};
    TowerInfoContainer *m_CaloWaveformContainer{nullptr};

    std::vector<std::vector<float>> m_waveforms = {{}};
    std::vector<std::vector<int>> *m_pedestal = nullptr;
    int m_runNumber{0};
    int m_nsamples{31};
    int m_nchannels{24576};
    float m_sampletime{50. / 3.};
    int m_fixpedestal{1500};
    int m_gaussian_noise{3};
    float m_deltaT{100.};
    float m_timeshiftwidth{0.};
    float m_peakpos{6.};
    gsl_rng *m_RandomGenerator{nullptr};

    PHG4CylinderCellGeom_Spacalv1 *geo{nullptr};
    const PHG4CylinderGeom_Spacalv3 *layergeom{nullptr};
    float m_sampling_fraction = {1.0};
    void maphitetaphi(PHG4Hit *g4hit, unsigned short &etabin, unsigned short &phibin, float &correction);
    unsigned int (*encode_tower)(const unsigned int etabin, const unsigned int phibin){TowerInfoDefs::encode_emcal};
    unsigned int (*decode_tower)(const unsigned int tower_key){TowerInfoDefs::decode_emcal};
    double template_function(double *x, double *par);
    void CreateNodeTree(PHCompositeNode *topNode);

    LightCollectionModel light_collection_model;

    NoiseType m_noiseType{NOISE_TREE};


};

#endif // CALOWAVEFORMSIM_H
