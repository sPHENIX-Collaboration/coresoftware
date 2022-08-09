#include "AnalysisExample.h"
#include <limits>
namespace HepMC3
{
HEPMC3_DECLARE_WRITER_FILE(AnalysisExample)
HEPMC3_DECLARE_WRITER_STREAM(AnalysisExample)

AnalysisExample::AnalysisExample(const std::string &filename,std::shared_ptr<GenRunInfo> /*run*/): m_file(filename),
    m_stream(&m_file)
{
    if ( !m_file.is_open() ) {
        HEPMC3_ERROR( "AnalysisExample: could not open output file: "<<filename )
    }
    m_sum_of_weights=0;
    m_sum_of_weights2=0;
    m_bins["rapidity"]=std::vector<double> {-std::numeric_limits<double>::infinity(), -5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,std::numeric_limits<double>::infinity()};
    m_vals["rapidity"]=std::vector<double>(m_bins.at("rapidity").size()-1,0.0);
    m_errs["rapidity"]=std::vector<double>(m_bins.at("rapidity").size()-1,0.0);
}

AnalysisExample::AnalysisExample(std::ostream &stream, std::shared_ptr<GenRunInfo> run)
    : m_file(),
      m_stream(&stream)
{
    m_sum_of_weights=0;
    m_sum_of_weights2=0;
    m_bins["rapidity"]=std::vector<double> {-std::numeric_limits<double>::infinity(), -5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,std::numeric_limits<double>::infinity()};
    m_vals["rapidity"]=std::vector<double>(m_bins.at("rapidity").size()-1,0.0);
    m_errs["rapidity"]=std::vector<double>(m_bins.at("rapidity").size()-1,0.0);

}

void AnalysisExample::write_event(const GenEvent &evt)
{
    double w=evt.weight();
    m_sum_of_weights+=w;
    m_sum_of_weights2+=w*w;
    for(auto p: evt.particles() )
    {
        if (p->status()!=1) continue;
        double eta=p->momentum().eta();
        int bin=std::distance(m_bins["rapidity"].begin(), lower_bound(m_bins["rapidity"].begin(),m_bins["rapidity"].end(),eta))-1;
        if (bin<0) bin=0;
        m_vals["rapidity"][bin]+=w;
        m_errs["rapidity"][bin]+=w*w;
    }
}

void AnalysisExample::close() {
    if (!m_stream) return;
    std::ofstream* ofs = dynamic_cast<std::ofstream*>(m_stream);
    for (size_t i=1; i<m_vals["rapidity"].size()-1; i++)
    {
        double val=m_vals["rapidity"][i]/m_sum_of_weights/(m_bins["rapidity"][i+1]-m_bins["rapidity"][i]);
        double err=sqrt(m_errs["rapidity"][i])/m_sum_of_weights/(m_bins["rapidity"][i+1]-m_bins["rapidity"][i]);
        (*ofs)<< std::fixed  << std::setprecision( 6 )<<m_bins["rapidity"][i]<<" "<<m_bins["rapidity"][i+1]<<" "<<val<<" "<<err<<std::endl;
    }
    if (ofs && !ofs->is_open()) return;
    if (ofs) ofs->close();

}

} // namespace HepMC3
