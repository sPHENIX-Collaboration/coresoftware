# Translated into python using testPolarization.cc (garren@fnal.gov, Oct. 2010) as an example
#
# andrii.verbytskyi@mpp.mpg.org, Nov. 2018
#

from pyHepMC3TestUtils import update_path, python_label
import sys

sys.path = update_path()

import random, math
from pyHepMC3TestUtils import COMPARE_ASCII_FILES
from pyHepMC3 import HepMC3 as hm


def test_Polarization():
    xout1 = hm.WriterAscii(python_label() + "testPolarization1.dat")
    xout2 = hm.WriterAscii(python_label() + "testPolarization2.dat")
    xout4 = hm.WriterAsciiHepMC2(python_label() + "testPolarization4.dat")
    xout5 = hm.WriterAscii(python_label() + "testPolarization5.dat")

    # Build the graph, which will look like
    # Please note this is not physically meaningful event.
    #                       p7                   #
    # p1                   /                     #
    #   \v1__p3      p5---v4                     #
    #         \_v3_/       \                     #
    #         /    \        p8                   #
    #    v2__p4     \                            #
    #   /            p6                          #
    # p3                                         #
    #
    # define a flow pattern as  p1 . p3 . p6
    #                       and p2 . p4 . p5
    #

    # First create the event container, with Signal Process 20, event number 1
    #
    evt = hm.GenEvent(hm.Units.GEV, hm.Units.MM)
    evt.set_event_number(1)
    evt.add_attribute("signal_process_id", hm.IntAttribute(20))

    v1 = hm.GenVertex()
    evt.add_vertex(v1)
    p1 = hm.GenParticle(hm.FourVector(0, 0, 7000, 7000), 2212, 3)
    evt.add_particle(p1)
    p1.add_attribute("flow1", hm.IntAttribute(231))
    p1.add_attribute("flow1", hm.IntAttribute(231))
    p1.add_attribute("theta", hm.DoubleAttribute(random.random() * math.pi))
    p1.add_attribute("phi", hm.DoubleAttribute(random.random() * math.pi * 2))
    v1.add_particle_in(p1)

    v2 = hm.GenVertex()
    evt.add_vertex(v2)
    p2 = hm.GenParticle(hm.FourVector(0, 0, -7000, 7000), 2212, 3)
    evt.add_particle(p2)
    p2.add_attribute("flow1", hm.IntAttribute(243))
    p2.add_attribute("theta", hm.DoubleAttribute(random.random() * math.pi))
    p2.add_attribute("phi", hm.DoubleAttribute(random.random() * math.pi * 2))
    v2.add_particle_in(p2)

    p3 = hm.GenParticle(hm.FourVector(0.751, -1.569, 32.191, 32.238), 1, 3)
    evt.add_particle(p3)
    p3.add_attribute("flow1", hm.IntAttribute(231))
    p3.add_attribute("theta", hm.DoubleAttribute(random.random() * math.pi))
    p3.add_attribute("phi", hm.DoubleAttribute(random.random() * math.pi * 2))
    v1.add_particle_out(p3)
    p4 = hm.GenParticle(hm.FourVector(-3.047, -19.0, -54.629, 57.920), -2, 3)
    evt.add_particle(p4)
    p4.add_attribute("flow1", hm.IntAttribute(243))
    p4.add_attribute("theta", hm.DoubleAttribute(random.random() * math.pi))
    p4.add_attribute("phi", hm.DoubleAttribute(random.random() * math.pi * 2))
    v2.add_particle_out(p4)

    v3 = hm.GenVertex()
    evt.add_vertex(v3)
    v3.add_particle_in(p3)
    v3.add_particle_in(p4)
    p6 = hm.GenParticle(hm.FourVector(-3.813, 0.113, -1.833, 4.233), 22, 1)
    evt.add_particle(p6)
    p6.add_attribute("flow1", hm.IntAttribute(231))
    p6.add_attribute("theta", hm.DoubleAttribute(random.random() * math.pi))
    p6.add_attribute("phi", hm.DoubleAttribute(random.random() * math.pi * 2))
    v3.add_particle_out(p6)
    p5 = hm.GenParticle(hm.FourVector(1.517, -20.68, -20.605, 85.925), -24, 3)
    evt.add_particle(p5)
    p5.add_attribute("flow1", hm.IntAttribute(243))
    p5.add_attribute("theta", hm.DoubleAttribute(random.random() * math.pi))
    p5.add_attribute("phi", hm.DoubleAttribute(random.random() * math.pi * 2))
    v3.add_particle_out(p5)
    # create v4
    v4 = hm.GenVertex(hm.FourVector(0.12, -0.3, 0.05, 0.004))
    evt.add_vertex(v4)
    v4.add_particle_in(p5)
    p7 = hm.GenParticle(hm.FourVector(-2.445, 28.816, 6.082, 29.552), 1, 1)
    evt.add_particle(p7)
    v4.add_particle_out(p7)
    p8 = hm.GenParticle(hm.FourVector(3.962, -49.498, -26.687, 56.373), -2, 1)
    evt.add_particle(p8)
    v4.add_particle_out(p8)

    evt.add_attribute("signal_process_vertex", hm.IntAttribute(v3.id()))
    evt.set_beam_particles(p1,p2)
    # The event is complete, we now print it out
    hm.Print.content(evt)
    hm.Print.listing(evt, 8)
    print(hm.version())
    print(hm.Print.line(v4, True))
    a = 0
    print(evt.particles())
    print(len(evt.particles()))
    print(evt.particles()[0])
    for ip in evt.particles():
        print(hm.Print.line(ip, True))
    xout1.write_event(evt)

    # write event in old format
    xout4.write_event(evt)
    # make a copy and write it
    xout5.write_event(hm.GenEvent(evt))
    # try changing polarization
    p2.add_attribute("theta", hm.DoubleAttribute(random.random() * math.pi))
    p2.add_attribute("phi", hm.DoubleAttribute(random.random() * math.pi))
    xout2.write_event(evt)
    xout1.close()
    xout2.close()
    xout4.close()
    xout5.close()
    # now clean-up by deleteing all objects from memory
    #
    # deleting the event deletes all contained vertices, and all particles
    # contained in those vertices
    evt.clear()

    assert (
        COMPARE_ASCII_FILES(python_label() + "testPolarization1.dat", python_label() + "testPolarization5.dat") == 0
    ) and (COMPARE_ASCII_FILES(python_label() + "testPolarization1.dat", python_label() + "testPolarization2.dat") != 0)

    inputA1 = hm.ReaderAscii(python_label() + "testPolarization1.dat")
    if inputA1.failed():
        return 1
    while not inputA1.failed():
        evt = hm.GenEvent()
        inputA1.read_event(evt)
        if inputA1.failed():
            print("End of file reached. Exit.\n")
            break
        evt.clear()
    inputA1.close()
    inputA2 = hm.ReaderAscii(python_label() + "testPolarization2.dat")
    if inputA2.failed():
        return 2
    while not inputA2.failed():
        evt = hm.GenEvent()
        inputA2.read_event(evt)
        if inputA2.failed():
            print("End of file reached. Exit.\n")
            break
        evt.clear()
    inputA2.close()

    inputA4 = hm.ReaderAsciiHepMC2(python_label() + "testPolarization4.dat")
    if inputA4.failed():
        return 4
    while not inputA4.failed():
        evt = hm.GenEvent()
        inputA4.read_event(evt)
        if inputA4.failed():
            print("End of file reached. Exit.\n")
            break
        evt.clear()
    inputA4.close()

    inputA5 = hm.ReaderAscii(python_label() + "testPolarization5.dat")
    if inputA5.failed():
        return 4
    while not inputA5.failed():
        evt = hm.GenEvent()
        inputA5.read_event(evt)
        if inputA5.failed():
            print("End of file reached. Exit.\n")
            break
        evt.clear()
    inputA5.close()

    return 0


if __name__ == "__main__":
    result = 1
    try:
        result = test_Polarization()
    except:
        result = 1
    sys.exit(result)
