#include "QAKFParticleTrackPtAsymmetry.h"

#include <qautils/QAHistManagerDef.h> // for getHistoManager
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <KFParticle.h>

#include <TH1.h>
#include <TH2.h>

#include <cassert>
#include <iostream>

QAKFParticleTrackPtAsymmetry::QAKFParticleTrackPtAsymmetry(const std::string &histo_prefix,            //
                                                           double min_m,                               //
                                                           double max_m,                               //
                                                           const std::vector<double> &mother_eta_bins, //
                                                           const std::vector<double> &mother_phi_bins)
    : m_prefix(histo_prefix)
    , m_min_mass(min_m)
    , m_max_mass(max_m)
    , m_mother_eta_bins(mother_eta_bins)
    , m_mother_phi_bins(mother_phi_bins)
{
    if (m_mother_eta_bins.size() < 2)
    {
        std::cout << __PRETTY_FUNCTION__ << " Error: need at least 2 eta bin edges. Use default" << std::endl;
        m_mother_eta_bins = std::vector<double>{-2.0, -1.0, -0.5, 0, 0.5, 1.0, 2.0};
    }

    if (m_mother_phi_bins.size() < 2)
    {
        std::cout << __PRETTY_FUNCTION__ << " Error: need at least 2 phi bin edges. Use default" << std::endl;
        m_mother_phi_bins = std::vector<double>{-3.2, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.2};
    }
}

void QAKFParticleTrackPtAsymmetry::bookHistograms(Fun4AllHistoManager *hm)
{
    assert(hm);

    // basic histograms
    h_trackPtAsymmetry = new TH1F(TString(m_prefix) + "TrackPtAsymmetry", ";#Delta p_{T}(trk1, trk2)/#Sigma p_{T}(trk1, trk2);Entries", 100, -1.0, 1.0);
    hm->registerHisto(h_trackPtAsymmetry);

    h2_trackPtAsymmetry_vs_mass = new TH2F(TString(m_prefix) + "TrackPtAsymmetry_vs_Mass", ";#Delta p_{T}(trk1, trk2)/#Sigma p_{T}(trk1, trk2);Mass [GeV]", 100, -1.0, 1.0, 100, m_min_mass, m_max_mass);
    hm->registerHisto(h2_trackPtAsymmetry_vs_mass);

    // differential histograms in (eta, phi) bins
    for (size_t i = 0; i < m_mother_eta_bins.size() - 1; ++i)
    {
        std::vector<TH2 *> phi_bin_histos;
        for (size_t j = 0; j < m_mother_phi_bins.size() - 1; ++j)
        {
            std::string histo_name = TString::Format("%sTrackPtAsymmetry_eta%.2f_%.2f_phi%.2f_%.2f", //
                                                     m_prefix.c_str(),                          //
                                                     m_mother_eta_bins[i],                      //
                                                     m_mother_eta_bins[i + 1],                  //
                                                     m_mother_phi_bins[j],                      //
                                                     m_mother_phi_bins[j + 1])
                                         .Data();
            TH2 *h2 = new TH2F(histo_name.c_str(),
                               TString::Format(";#Delta p_{T}(trk1, trk2)/#Sigma p_{T}(trk1, trk2);Mass [GeV]"),
                               100,
                               -1.0,
                               1.0,
                               100,
                               m_min_mass,
                               m_max_mass);
            hm->registerHisto(h2);
            phi_bin_histos.push_back(h2);
        }
        h2_trackPtAsymmetry_eta_phi_bins.push_back(phi_bin_histos);
    }
}

void QAKFParticleTrackPtAsymmetry::analyzeTrackPtAsymmetry(SvtxTrackMap *m_trackMap, KFParticle *mother)
{
    if (!m_trackMap)
    {
        if (m_verbosity > 0)
        {
            std::cout << __PRETTY_FUNCTION__ << " Error: null SvtxTrackMap pointer" << std::endl;
        }
        return;
    }

    if (!mother)
    {
        if (m_verbosity > 0)
        {
            std::cout << __PRETTY_FUNCTION__ << " Error: null mother pointer" << std::endl;
        }
        return;
    }
    else
    {
        if (mother->DaughterIds().size() != 2)
        {
            if (m_verbosity > 0)
            {
                std::cout << __PRETTY_FUNCTION__ << " Error: number of decay products is " << mother->DaughterIds().size() << ", not a 2-body decay" << std::endl;
            }
            return;
        }
    }

    const std::vector<int> track_ids = mother->DaughterIds();
    SvtxTrack *trk1 = nullptr;
    SvtxTrack *trk2 = nullptr;
    auto it1 = m_trackMap->find(track_ids[0]);
    if (it1 != m_trackMap->end())
    {
        trk1 = it1->second;
    }
    auto it2 = m_trackMap->find(track_ids[1]);
    if (it2 != m_trackMap->end())
    {
        trk2 = it2->second;
    }

    if (!trk1 || !trk2)
    {
        if (m_verbosity > 0)
        {
            std::cout << __PRETTY_FUNCTION__ << " Error: cannot to find one or both daughter tracks in SvtxTrackMap" << std::endl;
        }
        return;
    }

    //! positive and negative charged tracks (do we need this?)
    SvtxTrack *pos_trk = nullptr;
    SvtxTrack *neg_trk = nullptr;
    if (trk1->get_charge() > 0 && trk2->get_charge() < 0)
    {
        pos_trk = trk1;
        neg_trk = trk2;
    }
    else if (trk1->get_charge() < 0 && trk2->get_charge() > 0)
    {
        pos_trk = trk2;
        neg_trk = trk1;
    }
    else
    {
        if (m_verbosity > 0)
        {
            std::cout << __PRETTY_FUNCTION__ << " Error: both daughter tracks have the same charge" << std::endl;
        }
        return;
    }

    float pt_pos = pos_trk->get_pt();
    float pt_neg = neg_trk->get_pt();
    float pt_sum = pt_pos + pt_neg;
    if (pt_sum == 0)
    {
        if (m_verbosity > 0)
        {
            std::cout << __PRETTY_FUNCTION__ << " Warning: sum of daughter track pt is zero" << std::endl;
        }
        return;
    }
    float pt_asymmetry = (pt_neg - pt_pos) / pt_sum;

    if (h_trackPtAsymmetry)
    {
        h_trackPtAsymmetry->Fill(pt_asymmetry);
    }
    if (h2_trackPtAsymmetry_vs_mass)
    {
        h2_trackPtAsymmetry_vs_mass->Fill(pt_asymmetry, mother->GetMass());
    }

    float mother_eta = 0;
    float mother_eta_err = 0;
    mother->GetEta(mother_eta, mother_eta_err);
    float mother_phi = 0;
    float mother_phi_err = 0;
    mother->GetPhi(mother_phi, mother_phi_err);
    for (size_t i = 0; i < m_mother_eta_bins.size() - 1; ++i)
    {
        if (mother_eta >= m_mother_eta_bins[i] && mother_eta < m_mother_eta_bins[i + 1])
        {
            for (size_t j = 0; j < m_mother_phi_bins.size() - 1; ++j)
            {
                if (mother_phi >= m_mother_phi_bins[j] && mother_phi < m_mother_phi_bins[j + 1])
                {
                    TH2 *h2 = h2_trackPtAsymmetry_eta_phi_bins[i][j];
                    if (h2)
                    {
                        h2->Fill(pt_asymmetry, mother->GetMass());
                    }
                }
            }
        }
    }
}