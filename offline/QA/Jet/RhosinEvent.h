// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QA_JET_RHOSINEVENT_H
#define QA_JET_RHOSINEVENT_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <utility>  // std::pair, std::make_pair
#include <vector>

class Fun4AllHistoManager;
class PHCompositeNode;
class TH1;


/// \class RhosinEvent
/// \brief SubsysReco module that plots tower rho information
/// \author: Tanner Mengel

class RhosinEvent : public SubsysReco
{

    public:

        RhosinEvent(const std::string &name = "RhosinEvent");
        ~RhosinEvent() override {};

        void add_mult_rho_node(const std::string &name)
        {
            m_do_mult_rho = true;
            m_mult_rho_node = name;
        }
        void add_area_rho_node(const std::string &name)
        {
            m_do_area_rho = true;
            m_area_rho_node = name;
        }

        // standard Fun4All functions
        int Init(PHCompositeNode *topNode) override;
        int process_event(PHCompositeNode *topNode) override;
        int End(PHCompositeNode *topNode) override;


    private:

        //! Input Node strings
        // std::string m_outputFileName{"RhosinEvent.root"};

        bool m_do_mult_rho{true};
        bool m_do_area_rho{true};
        std::string m_mult_rho_node{"TowerRho_MULT"};
        std::string m_area_rho_node{"TowerRho_AREA"};

        Fun4AllHistoManager* m_manager{nullptr};
        //histos
        TH1 * h1_mult_rho{nullptr};
        TH1 * h1_mult_rho_sigma{nullptr};
        TH1 * h1_area_rho{nullptr};
        TH1 * h1_area_rho_sigma{nullptr};



};

#endif // QA_JET_RHOSINEVENT_H
