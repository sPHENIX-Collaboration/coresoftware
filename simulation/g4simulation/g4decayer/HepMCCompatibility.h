// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
///
/// @file HepMCCompatibility.h
/// @brief Implementation of compatibility layer (in-memory conversion functions) between HepMC2 and HepMC3
///
#ifndef HEPMCCOMPATIBILITY_H
#define HEPMCCOMPATIBILITY_H

#include <HepMC3/GenVertex.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenEvent.h>

///Please note the HEPMC_HAS_CENTRALITY should be defined externaly
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenEvent.h>

/** Converts HepMC3::Genevent to HepMC::Genevent */
HepMC::GenEvent* ConvertHepMCGenEvent_3to2(const HepMC3::GenEvent& evt)
{
    HepMC::GenEvent* n=new HepMC::GenEvent();
    HepMC::Units::LengthUnit   lu;
    HepMC::Units::MomentumUnit mu;
    if (evt.length_unit()==HepMC3::Units::CM) lu=HepMC::Units::CM;
    if (evt.length_unit()==HepMC3::Units::MM) lu=HepMC::Units::MM;
    if (evt.momentum_unit()==HepMC3::Units::MEV) mu=HepMC::Units::MEV;
    if (evt.momentum_unit()==HepMC3::Units::GEV) mu=HepMC::Units::GEV;
    n->use_units(mu,lu);


    std::shared_ptr<HepMC3::IntAttribute> A_signal_process_vertex=evt.attribute<HepMC3::IntAttribute>("signal_process_vertex");
    int signal_process_vertex_id=0;
    if (A_signal_process_vertex) signal_process_vertex_id=A_signal_process_vertex->value();

    std::map<HepMC3::ConstGenVertexPtr,HepMC::GenVertex*> vertexmap3to2;
    for (auto v3: evt.vertices())
    {
        HepMC3::FourVector pos3=v3->position();
        HepMC::FourVector pos2=HepMC::FourVector(pos3.x(),pos3.y(),pos3.pz(),pos3.t());
        HepMC::GenVertex* v2= new HepMC::GenVertex(pos2);
        std::vector<double> vweights;
        if(v3->attribute_names().size())
        {
            std::shared_ptr<HepMC3::VectorDoubleAttribute> rsvec=v3->attribute<HepMC3::VectorDoubleAttribute>("weights");
            if (rsvec) { vweights =rsvec->value(); }
            else {
                for (size_t ii=0; ii<100; ii++)
                {
                    std::shared_ptr<HepMC3::DoubleAttribute> rs=v3->attribute<HepMC3::DoubleAttribute>("weight"+std::to_string((long long unsigned int)ii));
                    if (!rs) break;
                    vweights.push_back(rs->value());
                }
            }
        }
        if (vweights.empty())vweights.push_back(1.0);
        v2->weights()=vweights;
        v2->suggest_barcode(v3->id());
        n->add_vertex(v2);
        if (signal_process_vertex_id==v3->id())  n->set_signal_process_vertex(v2);
        vertexmap3to2[v3]=v2;
    }
    std::map<HepMC3::ConstGenParticlePtr,HepMC::GenParticle*> particlemap3to2;
    for (auto p3: evt.particles())
    {
        HepMC3::FourVector mom3=p3->momentum();
        HepMC::FourVector mom2=HepMC::FourVector(mom3.px(),mom3.py(),mom3.pz(),mom3.e());
        HepMC::GenParticle* p2= new HepMC::GenParticle(mom2,p3->pid(),p3->status());
        if (p3->is_generated_mass_set())p2->set_generated_mass(p3->generated_mass());
        /** Converts HepMC3::Genevent to HepMC::Genevent */
        p2->suggest_barcode(10000+p3->id());
        particlemap3to2[p3]=p2;

        auto v3production=p3->production_vertex();
        HepMC::GenVertex* v2production=nullptr;
        if (v3production)
        {
            auto v2=vertexmap3to2.find(v3production);
            if (v2!=vertexmap3to2.end())  v2production=vertexmap3to2[v3production];
        }

        auto v3end=p3->end_vertex();
        HepMC::GenVertex* v2end=nullptr;
        if (v3end)
        {
            auto v2=vertexmap3to2.find(v3end);
            if (v2!=vertexmap3to2.end())  v2end=vertexmap3to2[v3end];
        }

        if (v2production) v2production->add_particle_out(p2);
        if (v2end) v2end->add_particle_in(p2);

        std::shared_ptr<HepMC3::DoubleAttribute> p3_theta=p3->attribute<HepMC3::DoubleAttribute>("theta");
        std::shared_ptr<HepMC3::DoubleAttribute> p3_phi=p3->attribute<HepMC3::DoubleAttribute>("phi");
        if (p3_phi && p3_theta)p2->set_polarization(HepMC::Polarization(p3_theta->value(),p3_phi->value()));
        std::shared_ptr<HepMC3::VectorIntAttribute> flows = p3->attribute<HepMC3::VectorIntAttribute>("flows");
        if (flows)
        {
            std::vector<int> flowvalues=flows->value();
            for (size_t i=0; i<flowvalues.size(); i++) p2->set_flow(i+1,flowvalues.at(i));
        } else {
            std::shared_ptr<HepMC3::IntAttribute> p3_flow1=p3->attribute<HepMC3::IntAttribute>("flow1");
            std::shared_ptr<HepMC3::IntAttribute> p3_flow2=p3->attribute<HepMC3::IntAttribute>("flow2");
            std::shared_ptr<HepMC3::IntAttribute> p3_flow3=p3->attribute<HepMC3::IntAttribute>("flow3");
            if (p3_flow1) p2->set_flow(1,p3_flow1->value());
            if (p3_flow2) p2->set_flow(2,p3_flow2->value());
            if (p3_flow3) p2->set_flow(3,p3_flow2->value());
        }
    }
    std::vector<HepMC3::ConstGenParticlePtr> bms=evt.beams();
    if (bms.size()>=2)
    {
        if (particlemap3to2.find(bms[0])!=particlemap3to2.end() && particlemap3to2.find(bms[1])!=particlemap3to2.end() )
        {
            n->set_beam_particles(particlemap3to2[bms[0]],particlemap3to2[bms[1]]);
        }
    }
    std::shared_ptr<HepMC3::GenCrossSection> cs3 = evt.attribute<HepMC3::GenCrossSection>("GenCrossSection");
    if(cs3) {
        HepMC::GenCrossSection cs2;
        cs2.set_cross_section(cs3->xsec(),cs3->xsec_err());
        n->set_cross_section(cs2);
    }
    std::shared_ptr<HepMC3::GenPdfInfo> pdf3 = evt.attribute<HepMC3::GenPdfInfo>("GenPdfInfo");
    if(pdf3)
    {
        HepMC::PdfInfo pdf2(
            pdf3->parton_id[0],
            pdf3->parton_id[1],
            pdf3->x[0],
            pdf3->x[1],
            pdf3->scale,
            pdf3->xf[0],
            pdf3->xf[1],
            pdf3->pdf_id[0],
            pdf3->pdf_id[1]
        );
        n->set_pdf_info(pdf2);
    }

    std::shared_ptr<HepMC3::GenHeavyIon> hi3 = evt.attribute<HepMC3::GenHeavyIon>("GenHeavyIon");
    if(hi3)
    {
        HepMC::HeavyIon hi2;
        hi2.set_Ncoll_hard(hi3->Ncoll_hard);
        hi2.set_Npart_proj(hi3->Npart_proj);
        hi2.set_Npart_targ(hi3->Npart_targ);
        hi2.set_Ncoll(hi3->Ncoll);
        hi2.set_spectator_neutrons(hi3->spectator_neutrons);
        hi2.set_spectator_protons(hi3->spectator_protons);
        hi2.set_N_Nwounded_collisions(hi3->N_Nwounded_collisions);
        hi2.set_Nwounded_N_collisions(hi3->Nwounded_N_collisions);
        hi2.set_Nwounded_Nwounded_collisions(hi3->Nwounded_Nwounded_collisions);
        hi2.set_impact_parameter(hi3->impact_parameter);
        hi2.set_event_plane_angle(hi3->event_plane_angle);
        hi2.set_eccentricity(hi3->eccentricity);
        hi2.set_sigma_inel_NN(hi3->sigma_inel_NN);
#ifdef HEPMC_HAS_CENTRALITY
        hi2.set_centrality(hi3->centrality);
#endif
        n->set_heavy_ion(hi2);
    }
	
	/*  
	
	Not using weight here for now to convert HepMC3 to HepMC
    std::vector<double> wv=evt.weights();
	std::vector<std::string> wn=evt.weight_names();
    for (size_t i=0; i<wv.size(); i++) n->weights()[wn.at(i)]=wv.at(i);
	
	*/
    

	std::shared_ptr<HepMC3::DoubleAttribute> A_event_scale=evt.attribute<HepMC3::DoubleAttribute>("event_scale");
    std::shared_ptr<HepMC3::DoubleAttribute> A_alphaQED=evt.attribute<HepMC3::DoubleAttribute>("alphaQED");
    std::shared_ptr<HepMC3::DoubleAttribute> A_alphaQCD=evt.attribute<HepMC3::DoubleAttribute>("alphaQCD");
    std::shared_ptr<HepMC3::IntAttribute> A_signal_process_id=evt.attribute<HepMC3::IntAttribute>("signal_process_id");
    std::shared_ptr<HepMC3::IntAttribute> A_mpi=evt.attribute<HepMC3::IntAttribute>("mpi");


    double event_scale=A_event_scale?(A_event_scale->value()):0.0;
    double alphaQED=A_alphaQED?(A_alphaQED->value()):0.0;
    double alphaQCD=A_alphaQCD?(A_alphaQCD->value()):0.0;
    int signal_process_id=A_signal_process_id?(A_signal_process_id->value()):0;
    int mpi=A_mpi?(A_mpi->value()):0;

    std::vector<long> random_states;
    for (int i=0; i<100; i++)
    {
        std::shared_ptr<HepMC3::IntAttribute> rs=evt.attribute<HepMC3::IntAttribute>("random_states"+std::to_string((long long unsigned int)i));
        if (!rs) break;
        random_states.push_back(rs->value());
    }
    n->set_mpi( mpi );
    n->set_signal_process_id( signal_process_id );
    n->set_event_number( evt.event_number() );
    n->set_random_states( random_states );
    n->set_event_scale( event_scale );
    n->set_alphaQCD( alphaQCD );
    n->set_alphaQED( alphaQED );
    return n;
}
/** Converts HepMC::Genevent to HepMC3::Genevent */
HepMC3::GenEvent*  ConvertHepMCGenEvent_2to3(const HepMC::GenEvent& evt, std::shared_ptr<HepMC3::GenRunInfo> run )
{
    HepMC3::GenEvent* n=new HepMC3::GenEvent();
    HepMC3::Units::LengthUnit lu;
    HepMC3::Units::MomentumUnit mu;
    if (evt.length_unit()==HepMC::Units::CM) lu=HepMC3::Units::CM;
    if (evt.length_unit()==HepMC::Units::MM) lu=HepMC3::Units::MM;
    if (evt.momentum_unit()==HepMC::Units::MEV) mu=HepMC3::Units::MEV;
    if (evt.momentum_unit()==HepMC::Units::GEV) mu=HepMC3::Units::GEV;
    n->set_units(mu,lu);

    n->set_run_info(run);
    n->set_event_number(evt.event_number());
    std::map<HepMC::GenVertex*,HepMC3::GenVertexPtr> vertexmap2to3;
    for (auto v2=evt.vertices_begin(); v2!=evt.vertices_end(); v2++)
    {
        HepMC::FourVector pos2=(*v2)->position();
        HepMC3::FourVector pos3=HepMC3::FourVector(pos2.x(),pos2.y(),pos2.pz(),pos2.t());
        HepMC3::GenVertexPtr v3=std::make_shared<HepMC3::GenVertex>(pos3);
        n->add_vertex(v3);
        vertexmap2to3[*v2]=v3;
    }

    std::map<HepMC::GenParticle*,HepMC3::GenParticlePtr> particlemap2to3;
    for (auto p2=evt.particles_begin(); p2!=evt.particles_end(); p2++)
    {

        HepMC::FourVector mom2=(*p2)->momentum();
        HepMC3::FourVector mom3=HepMC3::FourVector(mom2.px(),mom2.py(),mom2.pz(),mom2.e());
        HepMC3::GenParticlePtr p3= std::make_shared<HepMC3::GenParticle>(mom3,(*p2)->pdg_id(),(*p2)->status());
        /// we set it always as there is no way to check if it is set
        p3->set_generated_mass((*p2)->generated_mass());
        particlemap2to3[*p2]=p3;

        auto v2production=(*p2)->production_vertex();
        HepMC3::GenVertexPtr v3production;
        if (v2production)
        {
            auto v3=vertexmap2to3.find(v2production);
            if (v3!=vertexmap2to3.end())  v3production=vertexmap2to3[v2production];
        }

        auto v2end=(*p2)->end_vertex();
        HepMC3::GenVertexPtr v3end;
        if (v2end)
        {
            auto v3=vertexmap2to3.find(v2end);
            if (v3!=vertexmap2to3.end())  v3end=vertexmap2to3[v2end];
        }

        if (v3production) v3production->add_particle_out(p3);
        if (v3end) v3end->add_particle_in(p3);

        std::shared_ptr<HepMC3::DoubleAttribute> p3_theta=std::make_shared<HepMC3::DoubleAttribute>((*p2)->polarization().theta());

        std::shared_ptr<HepMC3::DoubleAttribute> p3_phi=std::make_shared<HepMC3::DoubleAttribute>((*p2)->polarization().phi());
        if (std::abs(p3_phi->value())>std::numeric_limits<double>::epsilon() || std::abs(p3_theta->value())>std::numeric_limits<double>::epsilon())
        {
            p3->add_attribute("theta",p3_theta);
            p3->add_attribute("phi",p3_theta);
        }

        std::vector<int> flv;
        auto fl=(*p2)->flow();
        for (auto flel=fl.begin(); flel!=fl.end(); ++flel) flv.push_back((*flel).second);
        std::shared_ptr<HepMC3::VectorIntAttribute> flows = std::make_shared<HepMC3::VectorIntAttribute>(flv);
        if (flv.size())
            p3->add_attribute("flows",flows);
    }


    std::shared_ptr<HepMC3::GenCrossSection> cs3 = std::make_shared<HepMC3::GenCrossSection>();
    auto cs2=evt.cross_section();
    if (cs2)cs3->set_cross_section(cs2->cross_section(),cs2->cross_section_error());
    n->set_cross_section(cs3);


    std::shared_ptr<HepMC3::GenPdfInfo> pdf3 = std::make_shared<HepMC3::GenPdfInfo>();
    auto pdf2=evt.pdf_info();
    if(pdf2)
    {
        pdf3->parton_id[0]=pdf2->id1();
        pdf3->parton_id[1]=pdf2->id2();
        pdf3->x[0]=pdf2->x1();
        pdf3->x[1]=pdf2->x2();
        pdf3->scale=pdf2->scalePDF();
        pdf3->xf[0]=pdf2->pdf1();
        pdf3->xf[1]=pdf2->pdf2();
        pdf3->pdf_id[0]=pdf2->pdf_id1();
        pdf3->pdf_id[1]=pdf2->pdf_id2();
        n->set_pdf_info(pdf3);
    }

    std::shared_ptr<HepMC3::GenHeavyIon> hi3 = std::make_shared<HepMC3::GenHeavyIon>();
    auto hi2=evt.heavy_ion();
    if(hi2)
    {
        hi3->Ncoll_hard=hi2->Ncoll_hard();
        hi3->Npart_proj=hi2->Npart_proj();
        hi3->Npart_targ=hi2->Npart_targ();
        hi3->Ncoll=hi2->Ncoll();
        hi3->spectator_neutrons=hi2->spectator_neutrons();
        hi3->spectator_protons=hi2->spectator_protons();
        hi3->N_Nwounded_collisions=hi2->N_Nwounded_collisions();
        hi3->Nwounded_N_collisions=hi2->Nwounded_N_collisions();
        hi3->Nwounded_Nwounded_collisions=hi2->Nwounded_Nwounded_collisions();
        hi3->impact_parameter=hi2->impact_parameter();
        hi3->event_plane_angle=hi2->event_plane_angle();
        hi3->eccentricity=hi2->eccentricity();
        hi3->sigma_inel_NN=hi2->sigma_inel_NN();
#ifdef HEPMC_HAS_CENTRALITY
        hi3->centrality=hi2->centrality();
#endif
        n->set_heavy_ion(hi3);
    }
    /** Yes, the desing is not always perfect */
    std::stringstream ss;
    evt.weights().print(ss);
    std::string allweights=ss.str();
    std::vector<std::string> wnames;
    std::string token;
    std::size_t pos=0;
    for (;;)
    {
        std::size_t name_begin=allweights.find_first_not_of(" (\n",pos);
        if (name_begin==std::string::npos) break;
        std::size_t name_end=allweights.find_first_of(",",name_begin);
        wnames.push_back(allweights.substr(name_begin,name_end-name_begin));
        pos=allweights.find_first_of(")",name_end)+1;
    }
    n->run_info()->set_weight_names(wnames);
    n->weights()=std::vector<double>(wnames.size(),1.0);
    for (auto wn: wnames)
    {
        n->weight(wn)=evt.weights()[wn];
    }

    std::shared_ptr<HepMC3::DoubleAttribute> A_event_scale=std::make_shared<HepMC3::DoubleAttribute>(evt.event_scale());
    n->add_attribute("event_scale",A_event_scale);
    std::shared_ptr<HepMC3::DoubleAttribute> A_alphaQED=std::make_shared<HepMC3::DoubleAttribute>(evt.alphaQED());
    n->add_attribute("alphaQED",A_alphaQED);
    std::shared_ptr<HepMC3::DoubleAttribute> A_alphaQCD=std::make_shared<HepMC3::DoubleAttribute>(evt.alphaQCD());
    n->add_attribute("alphaQCD",A_alphaQCD);
    std::shared_ptr<HepMC3::IntAttribute> A_signal_process_id=std::make_shared<HepMC3::IntAttribute>(evt.signal_process_id());
    n->add_attribute("signal_process_id",A_signal_process_id);
    std::shared_ptr<HepMC3::IntAttribute> A_mpi=std::make_shared<HepMC3::IntAttribute>(evt.mpi());
    n->add_attribute("mpi",A_mpi);
    std::shared_ptr<HepMC3::VectorLongIntAttribute> A_random_states= std::make_shared<HepMC3::VectorLongIntAttribute>( evt.random_states());
    n->add_attribute("random_states",A_random_states);

    auto spv2=evt.signal_process_vertex();
    if (vertexmap2to3.find(spv2)!=vertexmap2to3.end())
    {
        int signal_process_vertex_id3=vertexmap2to3[spv2]->id();
        std::shared_ptr<HepMC3::IntAttribute> A_signal_process_vertex=std::make_shared<HepMC3::IntAttribute>(signal_process_vertex_id3);
        n->add_attribute("signal_process_vertex",A_signal_process_vertex);
    }
    return n;
}
#endif
