// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/ReaderRootTree.h"
#include "HepMC3/WriterRootTree.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3TestUtils.h"
#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>
#include <TH1D.h>
using namespace HepMC3;
const Int_t kMaxparticles = 2000;
const Int_t kMaxvertices = 2000;
#ifndef DOXYGEN_SHOULD_SKIP_THIS
class SomeAnalysis
{
public :
    TChain         *fChain;   //!pointer to the analyzed TTree or TChain
    // Declaration of leaf types
//GenEventData *hepmc3_event;
    Int_t           event_number;
    Int_t           momentum_unit;
    Int_t           length_unit;
    Int_t           particles_;
    Int_t           particles_pid[kMaxparticles];   //[particles_]
    Int_t           particles_status[kMaxparticles];   //[particles_]
    Bool_t          particles_is_mass_set[kMaxparticles];   //[particles_]
    Double_t        particles_mass[kMaxparticles];   //[particles_]
    Double_t        particles_momentum_m_v1[kMaxparticles];   //[particles_]
    Double_t        particles_momentum_m_v2[kMaxparticles];   //[particles_]
    Double_t        particles_momentum_m_v3[kMaxparticles];   //[particles_]
    Double_t        particles_momentum_m_v4[kMaxparticles];   //[particles_]
    Int_t           vertices_;
    Int_t           vertices_status[kMaxvertices];   //[vertices_]
    Double_t        vertices_position_m_v1[kMaxvertices];   //[vertices_]
    Double_t        vertices_position_m_v2[kMaxvertices];   //[vertices_]
    Double_t        vertices_position_m_v3[kMaxvertices];   //[vertices_]
    Double_t        vertices_position_m_v4[kMaxvertices];   //[vertices_]
    vector<double>  weights;
    Double_t        event_pos_m_v1;
    Double_t        event_pos_m_v2;
    Double_t        event_pos_m_v3;
    Double_t        event_pos_m_v4;
    vector<int>     links1;
    vector<int>     links2;
    vector<int>     attribute_id;
    vector<string>  attribute_name;
    vector<string>  attribute_string;

    // List of branches
    TBranch        *b_hepmc3_event_event_number;   //!
    TBranch        *b_hepmc3_event_momentum_unit;   //!
    TBranch        *b_hepmc3_event_length_unit;   //!
    TBranch        *b_hepmc3_event_particles_;   //!
    TBranch        *b_particles_pid;   //!
    TBranch        *b_particles_status;   //!
    TBranch        *b_particles_is_mass_set;   //!
    TBranch        *b_particles_mass;   //!
    TBranch        *b_particles_momentum_m_v1;   //!
    TBranch        *b_particles_momentum_m_v2;   //!
    TBranch        *b_particles_momentum_m_v3;   //!
    TBranch        *b_particles_momentum_m_v4;   //!
    TBranch        *b_hepmc3_event_vertices_;   //!
    TBranch        *b_vertices_status;   //!
    TBranch        *b_vertices_position_m_v1;   //!
    TBranch        *b_vertices_position_m_v2;   //!
    TBranch        *b_vertices_position_m_v3;   //!
    TBranch        *b_vertices_position_m_v4;   //!
    TBranch        *b_hepmc3_event_weights;   //!
    TBranch        *b_hepmc3_event_event_pos_m_v1;   //!
    TBranch        *b_hepmc3_event_event_pos_m_v2;   //!
    TBranch        *b_hepmc3_event_event_pos_m_v3;   //!
    TBranch        *b_hepmc3_event_event_pos_m_v4;   //!
    TBranch        *b_hepmc3_event_links1;   //!
    TBranch        *b_hepmc3_event_links2;   //!
    TBranch        *b_hepmc3_event_attribute_id;   //!
    TBranch        *b_hepmc3_event_attribute_name;   //!
    TBranch        *b_hepmc3_event_attribute_string;   //!


    void Init(TChain *tree)
    {
        if (!tree) return;
        fChain = tree;
        fChain->SetMakeClass(1);

        fChain->SetBranchAddress("event_number", &event_number, &b_hepmc3_event_event_number);
        fChain->SetBranchAddress("momentum_unit", &momentum_unit, &b_hepmc3_event_momentum_unit);
        fChain->SetBranchAddress("length_unit", &length_unit, &b_hepmc3_event_length_unit);
        fChain->SetBranchAddress("particles", &particles_, &b_hepmc3_event_particles_);
        fChain->SetBranchAddress("particles.pid", particles_pid, &b_particles_pid);
        fChain->SetBranchAddress("particles.status", particles_status, &b_particles_status);
        fChain->SetBranchAddress("particles.is_mass_set", particles_is_mass_set, &b_particles_is_mass_set);
        fChain->SetBranchAddress("particles.mass", particles_mass, &b_particles_mass);
        fChain->SetBranchAddress("particles.momentum.m_v1", particles_momentum_m_v1, &b_particles_momentum_m_v1);
        fChain->SetBranchAddress("particles.momentum.m_v2", particles_momentum_m_v2, &b_particles_momentum_m_v2);
        fChain->SetBranchAddress("particles.momentum.m_v3", particles_momentum_m_v3, &b_particles_momentum_m_v3);
        fChain->SetBranchAddress("particles.momentum.m_v4", particles_momentum_m_v4, &b_particles_momentum_m_v4);
        fChain->SetBranchAddress("vertices", &vertices_, &b_hepmc3_event_vertices_);
        fChain->SetBranchAddress("vertices.status", vertices_status, &b_vertices_status);
        fChain->SetBranchAddress("vertices.position.m_v1", vertices_position_m_v1, &b_vertices_position_m_v1);
        fChain->SetBranchAddress("vertices.position.m_v2", vertices_position_m_v2, &b_vertices_position_m_v2);
        fChain->SetBranchAddress("vertices.position.m_v3", vertices_position_m_v3, &b_vertices_position_m_v3);
        fChain->SetBranchAddress("vertices.position.m_v4", vertices_position_m_v4, &b_vertices_position_m_v4);
        fChain->SetBranchAddress("weights", &weights, &b_hepmc3_event_weights);
        fChain->SetBranchAddress("event_pos.m_v1", &event_pos_m_v1, &b_hepmc3_event_event_pos_m_v1);
        fChain->SetBranchAddress("event_pos.m_v2", &event_pos_m_v2, &b_hepmc3_event_event_pos_m_v2);
        fChain->SetBranchAddress("event_pos.m_v3", &event_pos_m_v3, &b_hepmc3_event_event_pos_m_v3);
        fChain->SetBranchAddress("event_pos.m_v4", &event_pos_m_v4, &b_hepmc3_event_event_pos_m_v4);
        fChain->SetBranchAddress("links1", &links1, &b_hepmc3_event_links1);
        fChain->SetBranchAddress("links2", &links2, &b_hepmc3_event_links2);
        fChain->SetBranchAddress("attribute_id", &attribute_id, &b_hepmc3_event_attribute_id);
        fChain->SetBranchAddress("attribute_name", &attribute_name, &b_hepmc3_event_attribute_name);
        fChain->SetBranchAddress("attribute_string", &attribute_string, &b_hepmc3_event_attribute_string);

        fChain->SetBranchStatus("*",0);
        fChain->SetBranchStatus("particles",1);
        fChain->SetBranchStatus("particles.pid",1);
        fChain->SetBranchStatus("particles.status",1);
        fChain->SetBranchStatus("particles.momentum.m_v1",1);
        fChain->SetBranchStatus("particles.momentum.m_v2",1);

    }
    SomeAnalysis(const std::string& file)
    {
        TChain* TempChain= new TChain("hepmc3_tree");
        TempChain->Add(file.c_str());
        Init(TempChain);
    }
};
#endif

int main()
{
//Plain tree
    TH1D H1("H1","Pt of pions;Events/100MeV;P_{T},GeV",1000,0,100);
    SomeAnalysis* A= new SomeAnalysis("inputIO4.root");
    if (!A->fChain->GetEntries()) return 10001;
    for (int entry=0; entry<A->fChain->GetEntries(); entry++)
    {
        A->fChain->GetEntry(entry);
        for (int i=0; i<A->particles_; i++)
            if (A->particles_status[i]==1&&(std::abs(A->particles_pid[i])==211||std::abs(A->particles_pid[i])==11))
                H1.Fill(std::sqrt(A->particles_momentum_m_v1[i]*A->particles_momentum_m_v1[i]+A->particles_momentum_m_v2[i]*A->particles_momentum_m_v2[i]) );
    }
    delete A;
//GenEvent
    TH1D H2("H2","Pt of pions;Events/100MeV;P_{T},GeV",1000,0,100);
    ReaderRootTree inputA("inputIO4.root");
    if(inputA.failed()) return 10002;
    while( !inputA.failed() )
    {
        GenEvent evt(Units::GEV,Units::MM);
        inputA.read_event(evt);
        if( inputA.failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        for (ConstGenParticlePtr p: evt.particles())
            if ( std::abs(p->status()) == 1 && (std::abs(p->pdg_id()) == 211||std::abs(p->pdg_id()) == 11) )
                H2.Fill( p->momentum().perp());
        evt.clear();
    }
    inputA.close();
//Comparison
    int diff=0;
    for (int i=0; i<H1.GetNbinsX(); i++)
    {
        double eps=std::abs(H1.GetBinContent(i)-H2.GetBinContent(i));
        if (eps<1e-5) continue;
        std::cout<<"Bin: "<<i<<" "<<H1.GetBinContent(i)<<" "<<H2.GetBinContent(i)<<std::endl;
        diff++;
    }
    H1.Print("All");
    H2.Print("All");
    return diff;
}
