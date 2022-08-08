//////////////////////////////////////////////////////////////////////////
// testWeights.cc
//
// garren@fnal.gov, January 2010
// test Weights
//////////////////////////////////////////////////////////////////////////

#include <assert.h>
#include <iostream>
#include <string>
#include <vector>

#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/Print.h"
#include <stdexcept>
#include <limits>
using namespace HepMC3;
int main()
{
    GenEvent evt;
    std::shared_ptr<GenRunInfo> run = std::make_shared<GenRunInfo>();;
    evt.set_run_info(run);
    // original functionality
    evt.weights().push_back(2.0);
    evt.weights().push_back(4.56);
    assert( std::abs(evt.weights()[0] - 2.0) < std::numeric_limits<double>::epsilon() );
    assert( std::abs(evt.weights()[1] - 4.56) < std::numeric_limits<double>::epsilon() );
    assert( evt.weights().size() == 2 );
    assert( !evt.weights().empty() );

    std::vector<double> vec;
    for( int i = 0; i < 15; ++i )
    {
        double x = (double)i + 0.14*(double)i;
        vec.push_back( x );
    }
    evt.weights() = vec;
    assert( evt.weights().size() == 15 );
    evt.weights().pop_back();
    assert( evt.weights().size() == 14 );

    // new functionality
    std::vector<std::string> names;
    for( size_t i = 0; i < evt.weights().size() - 1; ++i ) names.push_back(std::to_string((unsigned long long)i));
    std::string nm = "tau";
    names.push_back(nm);
    run->set_weight_names(names);

    evt.weight(nm) = 3.1;
    //assert( evt.weights().size() == (vs) );

    // lookup a nonexistent name
    try
    {
        double x = evt.weight("bad");
        std::cout << "lookup of nonexistent name returns " << x << std::endl;
    }
    catch (std::exception& e)
    {
        std::cout << e.what() << std::endl;
        std::cout << "HepMC testWeights: the above error is intentional" << std::endl;
    }
    Print::listing(evt);
    return 0;
}
