// -*- C++ -*-
//
#include "WriterRootTreeOPAL.h"
#include "TTree.h"
namespace HepMC3
{
WriterRootTreeOPAL::WriterRootTreeOPAL(const std::string &filename,std::shared_ptr<GenRunInfo> run):WriterRootTree::WriterRootTree(filename,"h10","h10",run) {}
void WriterRootTreeOPAL::init_branches()
{
    m_tree->Branch("Irun", &m_Irun);
    m_tree->Branch("Ievnt", &m_Ievnt);
    m_tree->Branch("Ebeam",&m_Ebeam);
}
void WriterRootTreeOPAL::write_event(const GenEvent &evt)
{
    m_Ievnt=evt.event_number();
    std::vector<size_t> beams;
    for (size_t i=0; i<evt.particles().size(); i++)
        if (evt.particles().at(i)->status()==4&&std::abs(evt.particles().at(i)->pid())==11)
            beams.push_back(i);

    if (beams.size()==2)
        m_Ebeam=std::abs(evt.particles().at(beams[0])->momentum().e());
    else
        m_Ebeam=std::abs(evt.particles().at(0)->momentum().e());
    WriterRootTree::write_event(evt);
}
void WriterRootTreeOPAL::set_run_number(const int nr) {m_Irun=nr;}
} // namespace HepMC3
