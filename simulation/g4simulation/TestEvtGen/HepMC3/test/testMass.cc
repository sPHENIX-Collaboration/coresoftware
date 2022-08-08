//-------------------------------------------------------------------
// testMass.cc.in
//
// garren@fnal.gov, March 2006
// Read events written by example_MyPythia.cc
// Select events containing a photon of pT > 25 GeV
// Add arbitrary PDF information to one of the good events
// Add arbitrary HeavyIon information to one of the good events
// Write the selected events and read them back in using an istream
//-------------------------------------------------------------------

#include <cmath>	// for min()
#include <ostream>

#include "HepMC3/GenParticle.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenPdfInfo.h"
#include "HepMC3/GenHeavyIon.h"


#include "HepMC3/Version.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/WriterAsciiHepMC2.h"

// define methods and classes used by this test
#include "IsGoodEvent.h"
using namespace HepMC3;
bool   massInfo( const GenEvent&, std::ostream& );

int main()
{
    // declare an input strategy to read the data produced with the
    ReaderAsciiHepMC2 ascii_in("inputMass.hepmc");
    if (ascii_in.failed()) return 1;
    // declare another IO_GenEvent for output
    WriterAsciiHepMC2 ascii_out("testMass1.out");
    // declare an instance of the event selection predicate
    IsGoodEventDIS is_good_event;
    // send version to output
    HepMC3::version();
    //........................................EVENT LOOP
    int icount=0;
    int num_good_events=0;
    double x1=0., x2=0., q=0., xf1=0., xf2=0.;
    GenEvent evt;
    while ( !ascii_in.failed() )
    {
        bool readOK=ascii_in.read_event(evt);
        if (!readOK) return 1;
        icount++;
        if ( icount%50==1 ) std::cout << "Processing Event Number " << icount<< " its # " << evt.event_number() << std::endl;
        if ( is_good_event(evt) )
        {
            if (num_good_events == 0 )
            {
                // add some arbitrary PDF information
                x1 = std::min(0.8, 0.07 * icount);
                x2 = 1-x1;
                q = 1.69 * icount;
                // use beam momentum
                if( evt.beams().size()==2 )
                {
                    GenParticlePtr bp1 = evt.beams().at(0);
                    GenParticlePtr bp2 = evt.beams().at(1);
                    xf1 = x1*bp1->momentum().p3mod();
                    xf2 = x2*bp1->momentum().p3mod();
                }
                else
                {
                    xf1 = x1*0.34;
                    xf2 = x2*0.34;
                }
                // provide optional pdf set id numbers (two ints at the end of the constructor)
                std::shared_ptr< GenPdfInfo> pdf = std::make_shared< GenPdfInfo>();
                evt.add_attribute("GenPdfInfo",pdf);
                pdf->set( 2, 3, x1, x2, q, xf1, xf2, 230, 230);
                // add some arbitrary HeavyIon information
                std::shared_ptr< GenHeavyIon> ion = std::make_shared< GenHeavyIon>();
                evt.add_attribute("GenHeavyIon",ion);
                ion->set(23,11,12,15,3,5,0,0,0,0.0145,0.0,0.0,0.0,0.23);
            }
            std::cout << "saving Event " << evt.event_number() << std::endl;
            if( evt.weights().size() > 0 )
            {
                std::cout << "Weights: ";
                for ( std::vector<double>::const_iterator w=evt.weights().begin(); w!=evt.weights().end(); ++w )
                    std::cout <<" "<<*w;
                std::cout << std::endl;
            }
            ascii_out.write_event(evt);
            ++num_good_events;
        }
        // clean up and get next event
        evt.clear();
    }
    //........................................PRINT RESULT
    std::cout << num_good_events << " out of " << icount
              << " processed events passed the cuts. Finished." << std::endl;
    ascii_in.close();
    ascii_out.close();
    // now read the file we just created
    // declare an input strategy
    std::ifstream istr( "testMass1.out" );
    if( !istr )
    {
        std::cerr << "testMass: cannot open " << std::endl;
        return 1;
    }
    ReaderAsciiHepMC2 xin(istr);
    if (xin.failed()) return 1;
    // declare another IO_GenEvent for output
    WriterAsciiHepMC2 xout("testMass2.out");
    if (xout.failed()) return 1;
    //........................................EVENT LOOP
    int ixin=0;
    while ( !xin.failed() )
    {
        bool readOK=xin.read_event(evt);
        if (!readOK) return 1;
        ixin++;
        std::cout << "reading Event " << evt.event_number() << std::endl;
        if( evt.weights().size() > 0 )
        {
            std::cout << "Weights: ";
            for ( std::vector<double>::const_iterator w=evt.weights().begin(); w!=evt.weights().end(); ++w )
                std::cout <<" "<<*w;
            std::cout << std::endl;
        }
        xout.write_event(evt);
        // look at mass info
        if (!  massInfo(evt,std::cout)) return 1;
        // clean up and get next event
        evt.clear();
    }
    //........................................PRINT RESULT
    std::cout << ixin << " events in the second pass. Finished." << std::endl;
    xin.close();
    xout.close();
    return 0;
}

bool massInfo( const GenEvent& e, std::ostream& os )
{
    for (ConstGenParticlePtr  p: e.particles()) {
        double gm = p->generated_mass();
        double m = p->momentum().m();
        double d = std::abs(m-gm);
        if( d > 1.0e-4 && gm>1.0e-4)
        {
            os << "Event " << e.event_number()
               << " Particle " << (p)->pdg_id()
               << " generated mass " << gm
               << " mass from momentum " << m
               << " difference " << d << std::endl;
            return false;
        }
    }
    return true;
}
