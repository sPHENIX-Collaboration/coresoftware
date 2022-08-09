from pyHepMC3 import HepMC3 as hm
import sys


class Pythia8ToHepMC3:
    def __init__(self):
        self.m_internal_event_number = 0
        self.m_print_inconsistency = True
        self.m_free_parton_warnings = True
        self.m_crash_on_problem = False
        self.m_convert_gluon_to_0 = False
        self.m_store_pdf = True
        self.m_store_proc = True
        self.m_store_xsec = True
        self.m_store_weights = True

    # The recommended method to convert Pythia events into HepMC ones
    def fill_next_event1(self, pythia, evt, ievnum):
        return self.fill_next_event(pythia.event, evt, ievnum, pythia.infoPython(), pythia.settings)

    # Alternative method to convert Pythia events into HepMC ones
    def fill_next_event(self, pyev, evt, ievnum, pyinfo, pyset):
        # 1. Error if no event passed.
        if evt is None:
            print("Pythia8ToHepMC3::fill_next_event error - passed null event.")
            return False
        # Event number counter.
        if ievnum >= 0:
            evt.set_event_number(ievnum)
            self.m_internal_event_number = ievnum
        else:
            evt.set_event_number(self.m_internal_event_number)
            self.m_internal_event_number = self.m_internal_event_number + 1
        evt.set_units(hm.Units.GEV, hm.Units.MM)
        #        // 2. Fill particle information
        hepevt_particles = []
        for i in range(0, pyev.size()):
            hepevt_particles.append(
                hm.GenParticle(
                    hm.FourVector(pyev[i].px(), pyev[i].py(), pyev[i].pz(), pyev[i].e()),
                    pyev[i].id(),
                    pyev[i].statusHepMC(),
                )
            )
            hepevt_particles[i].set_generated_mass(pyev[i].m())
        #        // 3. Fill vertex information and find beam particles.
        # For type compatibility
        vertex_cache = hm.GenEvent().vertices()
        beam_particles = hm.GenEvent().particles()
        for i in range(0, pyev.size()):
            mothers = pyev[i].motherList()
            if len(mothers) != 0:
                prod_vtx = hepevt_particles[mothers[0]].end_vertex()
                if prod_vtx is None:
                    prod_vtx = hm.GenVertex()
                    vertex_cache.append(prod_vtx)
                    for j in range(0, len(mothers)):
                        prod_vtx.add_particle_in(hepevt_particles[mothers[j]])
                prod_pos = hm.FourVector(pyev[i].xProd(), pyev[i].yProd(), pyev[i].zProd(), pyev[i].tProd())
                #                // Update vertex position if necessary
                if (not prod_pos.is_zero()) and prod_vtx.position().is_zero():
                    prod_vtx.set_position(prod_pos)
                prod_vtx.add_particle_out(hepevt_particles[i])
            else:
                beam_particles.append(hepevt_particles[i])
        #        // Add particles and vertices in topological order
        if len(beam_particles) < 2:
            print("There are  ", len(beam_particles), "!=2 particles without mothers")
            if self.m_crash_on_problem:
                sys.exit(1)
        evt.add_tree(beam_particles)
        #        //Attributes should be set after adding the particles to event
        for i in range(0, pyev.size()):
            #            // Colour flow uses index 1 and 2.
            colType = pyev[i].colType()
            if colType == -1 or colType == 1 or colType == 2:
                flow1 = 0
                flow2 = 0
                if colType == 1 or colType == 2:
                    flow1 = pyev[i].col()
                if colType == -1 or colType == 2:
                    flow2 = pyev[i].acol()
                hepevt_particles[i].add_attribute("flow1", hm.IntAttribute(flow1))
                hepevt_particles[i].add_attribute("flow2", hm.IntAttribute(flow2))
        #        // If hadronization switched on then no final coloured particles.
        if pyset == None:
            doHadr = self.m_free_parton_warnings and pyset.flag("HadronLevel:Hadronize")
        else:
            doHadr = pyset.flag("HadronLevel:all") and pyset.flag("HadronLevel:Hadronize")

        #        // 4. Check for particles which come from nowhere, i.e. are without
        #        // mothers or daughters. These need to be attached to a vertex, or else
        #        // they will never become part of the event.
        for i in range(1, pyev.size()):

            #            // Check for particles not added to the event
            #            // NOTE: We have to check if this step makes any sense in HepMC event standard
            if not hepevt_particles[i]:
                print("hanging particle ", i)
                prod_vtx = hm.GenVertex()
                prod_vtx.add_particle_out(hepevt_particles[i])
                evt.add_vertex(prod_vtx)

            #            // Also check for free partons (= gluons and quarks; not diquarks?).
            if doHadr and self.m_free_parton_warnings:
                if hepevt_particles[i].pid() == 21 and (hepevt_particles[i].end_vertex() is None):
                    print("gluon without end vertex ", i)
                    if self.m_crash_on_problem:
                        sys.exit(1)
                if abs(hepevt_particles[i].pid()) <= 6 and (hepevt_particles[i].end_vertex() is None):
                    print("quark without end vertex ", i)
                    if self.m_crash_on_problem:
                        sys.exit(1)

        #        // 5. Store PDF, weight, cross section and other event information.
        #        // Flavours of incoming partons.
        if self.m_store_pdf and pyinfo is not None:
            id1pdf = pyinfo.id1pdf()
            id2pdf = pyinfo.id2pdf()
            if self.m_convert_gluon_to_0:
                if id1pdf == 21:
                    id1pdf = 0
                if id2pdf == 21:
                    id2pdf = 0
            pdfinfo = hm.GenPdfInfo()
            pdfinfo.set(id1pdf, id2pdf, pyinfo.x1pdf(), pyinfo.x2pdf(), pyinfo.QFac(), pyinfo.pdf1(), pyinfo.pdf2())
            #            // Store PDF information.
            evt.set_pdf_info(pdfinfo)

        #        // Store process code, scale, alpha_em, alpha_s.
        if self.m_store_proc and pyinfo is not None:
            evt.add_attribute("mpi", hm.IntAttribute(pyinfo.nMPI()))
            evt.add_attribute("signal_process_id", hm.IntAttribute(pyinfo.code()))
            evt.add_attribute("event_scale", hm.DoubleAttribute(pyinfo.QRen()))
            evt.add_attribute("alphaQCD", hm.DoubleAttribute(pyinfo.alphaS()))
            evt.add_attribute("alphaQED", hm.DoubleAttribute(pyinfo.alphaEM()))

        #        // Store cross-section information in pb.
        if self.m_store_xsec and pyinfo is not None:
            xsec = hm.GenCrossSection()
            xsec.set_cross_section(pyinfo.sigmaGen() * 1e9, pyinfo.sigmaErr() * 1e9)
            evt.set_cross_section(xsec)

        #        // Store event weights.
        if self.m_store_weights and pyinfo is not None:
            evt.weights().clear()
            for iweight in range(0, pyinfo.nWeights()):
                evt.weights().append(pyinfo.weight(iweight))

        #        // Done.
        return True
