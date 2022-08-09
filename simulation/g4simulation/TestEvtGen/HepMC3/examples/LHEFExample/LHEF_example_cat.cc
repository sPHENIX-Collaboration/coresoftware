/**
 *  @example LHEF_example_cat.cc
 *  @brief Basic example of use of LHEF for reading and writing LHE files
 */
#include "HepMC3/LHEFAttributes.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/AssociatedParticle.h"
#include <iomanip>

using namespace HepMC3;

int main(int /*argc*/, char ** /*argv*/) {

    // Create Reader to read the example LHE file.
    LHEF::Reader reader("LHEF_example.lhe");

    // Create a HEPRUP attribute and initialize it from the reader.
    std::shared_ptr<HEPRUPAttribute> hepr = std::make_shared<HEPRUPAttribute>();
    hepr->heprup = reader.heprup;

    // There may be some XML tags in the LHE file which are
    // non-standard, but we can save them as well.
    hepr->tags = LHEF::XMLTag::findXMLTags(reader.headerBlock + reader.initComments);

    // Nowwe want to create a GenRunInfo object for the HepMC file, and
    // we add the LHEF attribute to that.
    std::shared_ptr<GenRunInfo> runinfo = std::make_shared<GenRunInfo>();
    runinfo->add_attribute("HEPRUP", hepr);

    // This is just a test to make sure we can add other attributes as
    // well.
    runinfo->add_attribute("NPRUP",
                           std::make_shared<FloatAttribute>(hepr->heprup.NPRUP));

    // We want to be able to convey the different event weights to
    // HepMC. In particular we need to add the names of the weights to
    // the GenRunInfo object.
    std::vector<std::string> weightnames;
    weightnames.push_back("0"); // The first weight is always the
                                // default weight with name "0".
    for ( int i = 0, N = hepr->heprup.weightinfo.size(); i < N; ++i )
        weightnames.push_back(hepr->heprup.weightNameHepMC(i));
    runinfo->set_weight_names(weightnames);

    // We also want to convey the information about which generators was
    // used to HepMC.
    for ( int i = 0, N = hepr->heprup.generators.size(); i < N; ++i ) {
        GenRunInfo::ToolInfo tool;
        tool.name =  hepr->heprup.generators[i].name;
        tool.version =  hepr->heprup.generators[i].version;
        tool.description =  hepr->heprup.generators[i].contents;
        runinfo->tools().push_back(tool);
    }

    // Now we want to start reading events from the LHE file and
    // translate them to HepMC.
    WriterAscii output("LHEF_example.hepmc3", runinfo);
    int neve = 0;
    while ( reader.readEvent() ) {
        ++neve;

        // To each GenEvent we want to add an attribute corresponding to
        // the HEPEUP. Also here there may be additional non-standard
        // information outside the LHEF <event> tags, which we may want to
        // add.
        std::shared_ptr<HEPEUPAttribute> hepe = std::make_shared<HEPEUPAttribute>();
        if ( reader.outsideBlock.length() )
            hepe->tags = LHEF:: XMLTag::findXMLTags(reader.outsideBlock);
        hepe->hepeup = reader.hepeup;
        GenEvent ev(runinfo, Units::GEV, Units::MM);
        ev.set_event_number(neve);

        // This is just a text to check that we can add additional
        // attributes to each event.
        ev.add_attribute("HEPEUP", hepe);
        ev.add_attribute("AlphaQCD",
                        std:: make_shared<DoubleAttribute>(hepe->hepeup.AQCDUP));
        ev.add_attribute("AlphaEM",
                         std::make_shared<DoubleAttribute>(hepe->hepeup.AQEDUP));
        ev.add_attribute("NUP",
                         std::make_shared<IntAttribute>(hepe->hepeup.NUP));
        ev.add_attribute("IDPRUP",
                         std::make_shared<LongAttribute>(hepe->hepeup.IDPRUP));

        // Now add the Particles from the LHE event to HepMC
        std::vector<GenParticlePtr> particles;
        std::map< std::pair<int,int>, GenVertexPtr> vertices;
        for ( int i = 0; i < hepe->hepeup.NUP; ++i )
        {
            particles.push_back(std::make_shared<GenParticle>(hepe->momentum(i),hepe->hepeup.IDUP[i],hepe->hepeup.ISTUP[i]));
            if (i<2) continue;
            std::pair<int,int> vertex_index(hepe->hepeup.MOTHUP[i].first,hepe->hepeup.MOTHUP[i].second);
            if (vertices.find(vertex_index)==vertices.end())vertices[vertex_index]=std::make_shared<GenVertex>();
            vertices[vertex_index]->add_particle_out(particles.back());
        }
        for ( auto v: vertices )
        {
            std::pair<int,int> vertex_index=v.first;
            GenVertexPtr          vertex=v.second;
            for (int i=vertex_index.first-1; i<vertex_index.second; i++) if (i>=0&&i<(int)particles.size()) vertex->add_particle_in(particles[i]);
        }
        for ( auto v: vertices ) ev.add_vertex(v.second);

        // And we also want to add the weights.
        std::vector<double> wts;
        for ( int i = 0, N = hepe->hepeup.weights.size(); i < N; ++i )
            wts.push_back(hepe->hepeup.weights[i].first);
        ev.weights() = wts;

        // Let's see if we can associate p1 and p2.
        ev.add_attribute("OtherIncoming",
                         std::make_shared<AssociatedParticle>(particles[1]), particles[0]->id());


        // And then we are done and can write out the GenEvent.
        output.write_event(ev);

    }

    output.close();

    // Now we wnat to make sure we can read in the HepMC file and
    // recreate the same info. To check this we try to reconstruct the
    // LHC file we read in above.
    ReaderAscii input("LHEF_example.hepmc3");
    LHEF::Writer writer("LHEF_example_out.lhe");
    hepr = std::shared_ptr<HEPRUPAttribute>();

    // The loop over all events in the HepMC file.
    while ( true ) {

        // Read in the next event.
        GenEvent ev(Units::GEV, Units::MM);
        if ( !input.read_event(ev) || ev.event_number() == 0 ) break;

        // Check that the first incoming particle still refers to the second.
        std::shared_ptr<AssociatedParticle> assoc =
          ev.attribute<AssociatedParticle>("OtherIncoming", 1);
        if ( !assoc || !assoc->associated() ||
             assoc->associated() != ev.particles()[1] ) return 3;

        // Make sure the weight names are the same.
        if ( input.run_info()->weight_names() != weightnames ) return 2;

        // For the first event we also go in and reconstruct the HEPRUP
        // information, and write it out to the new LHE file.
        if ( !hepr ) {
            hepr = ev.attribute<HEPRUPAttribute>("HEPRUP");

            // Here we also keep track of the additional non-standard info
            // we found in the original LHE file.
            for ( int i = 0, N = hepr->tags.size(); i < N; ++i )
                if ( hepr->tags[i]->name != "init" )
                    hepr->tags[i]->print(writer.headerBlock());

            // This is just a test that we can access other attributes
            // included in the GenRunInfo.
            hepr->heprup.NPRUP =
                int(input.run_info()->
                    attribute<FloatAttribute>("NPRUP")->value());

            // Then we write out the HEPRUP object.
            writer.heprup = hepr->heprup;
            if ( writer.heprup.eventfiles.size() >= 2 ) {
              writer.heprup.eventfiles[0].filename = "LHEF_example_1_out.plhe";
              writer.heprup.eventfiles[1].filename = "LHEF_example_2_out.plhe";
            }
            writer.init();

        }

        // Now we can access the HEPEUP attribute of the current event.
        std::shared_ptr<HEPEUPAttribute> hepe =
            ev.attribute<HEPEUPAttribute>("HEPEUP");

        // Again, there may be addisional non-standard information we want
        // to keep.
        for ( int i = 0, N = hepe->tags.size(); i < N; ++i )
            if ( hepe->tags[i]->name != "event" &&
                 hepe->tags[i]->name != "eventgroup" )
                hepe->tags[i]->print(writer.eventComments());

        // This is just a test that we can access other attributes
        // included in the GenRunInfo.
        hepe->hepeup.AQCDUP =
            ev.attribute<DoubleAttribute>("AlphaQCD")->value();
        hepe->hepeup.AQEDUP =
            ev.attribute<DoubleAttribute>("AlphaEM")->value();
        hepe->hepeup.NUP =
            ev.attribute<IntAttribute>("NUP")->value();
        hepe->hepeup.IDPRUP =
            ev.attribute<LongAttribute>("IDPRUP")->value();

        // And then we can write out the HEPEUP object.
        writer.hepeup = hepe->hepeup;
        writer.hepeup.heprup =  &writer.heprup;
        writer.writeEvent();

    }

}
