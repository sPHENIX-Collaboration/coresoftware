//
// Thanks to Bob McElrath and Frank Siegert for this test
// andrii.verbytskyi@mpp.mpg.gov, Nov. 2018

#include <fstream>

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"
#include "HepMC3/PrintStreams.h"
using namespace HepMC3;
int main()
{
    GenEvent p_event(Units::GEV, Units::MM);
    for(int i=0; i<10; i++)
    {
        FourVector vector(1.0,1.0,1.0,1.0);
        GenVertexPtr vertex=std::make_shared<GenVertex>();
        vertex->set_position(vector);
        vertex->set_id(i);
        for(int j=0; j<3; j++)
        {
            GenParticlePtr particle = std::make_shared<GenParticle>(vector,1,2);
            vertex->add_particle_in(particle);
        }
        for(int j=0; j<3; j++)
        {
            GenParticlePtr particle = std::make_shared<GenParticle>(vector,1,2);
            vertex->add_particle_out(particle);
        }
        p_event.add_vertex(vertex);
    }
    Print::listing(p_event);
    //
    Print::content(p_event);
    std::cout<<p_event;
    // cleanup
    p_event.clear();
    return 0;
}
