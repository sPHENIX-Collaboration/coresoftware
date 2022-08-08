/**
 *  @example pythia8_example.cc
 *  @brief Basic example of use for pythia8 interface
 *
 */
#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include "Pythia8/Pythia.h"
#ifdef HEPMC3_USE_INTERFACE_FROM_PYTHIA8
#include "Pythia8Plugins/HepMC3.h"
#else
#include "Pythia8ToHepMC3.h"
#endif
#include <iostream>
using namespace HepMC3;

/** Main program */
int main(int argc, char **argv) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <pythia_config_file> <output_hepmc3_file>" << std::endl;
        exit(-1);
    }

    Pythia8::Pythia pythia;
    Pythia8ToHepMC3 pythiaToHepMC;
    pythia.readFile(argv[1]);
    pythia.init();
    std::shared_ptr<GenRunInfo> run = std::make_shared<GenRunInfo>();
    struct GenRunInfo::ToolInfo generator={std::string("Pythia8"),std::to_string(PYTHIA_VERSION).substr(0,5),std::string("Used generator")};
    run->tools().push_back(generator);
    struct GenRunInfo::ToolInfo config={std::string(argv[1]),"1.0",std::string("Control cards")};
    run->tools().push_back(config);
    std::vector<std::string> names;
    for (int iWeight=0; iWeight < pythia.info.nWeights(); ++iWeight) {
     std::string s=pythia.info.weightLabel(iWeight);
     if (!s.length()) s=std::to_string((long long int)iWeight);
     names.push_back(s);
    }
    if (!names.size()) names.push_back("default");
    run->set_weight_names(names);
    WriterAscii file(argv[2],run);

    int nEvent = pythia.mode("Main:numberOfEvents");

    for( int i = 0; i< nEvent; ++i ) {
        if( !pythia.next() ) continue;

        GenEvent hepmc( Units::GEV, Units::MM );

        pythiaToHepMC.fill_next_event(pythia.event, &hepmc, -1, &pythia.info);

        if( i==0 ) {
            std::cout << "First event: " << std::endl;
            Print::listing(hepmc);
        }

        file.write_event(hepmc);
    }

    file.close();
    pythia.stat();
}
