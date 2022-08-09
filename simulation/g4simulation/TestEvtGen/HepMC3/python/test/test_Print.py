from pyHepMC3TestUtils import update_path, python_label
import sys

sys.path = update_path()

import math, random
from pyHepMC3TestUtils import COMPARE_ASCII_FILES
from pyHepMC3 import HepMC3 as hm
import io


def test_Print():
    evt = hm.GenEvent(hm.Units.MomentumUnit.GEV, hm.Units.LengthUnit.CM)
    evt.set_event_number(1)
    evt.add_attribute("signal_process_id", hm.IntAttribute(20))
    #     create vertex 1
    v1 = hm.GenVertex()
    evt.add_vertex(v1)
    p1 = hm.GenParticle(hm.FourVector(1.0, 1.0, 7000, 7000), 2212, 3)
    evt.add_particle(p1)
    p1.add_attribute("flow1", hm.IntAttribute(231))
    p1.add_attribute("flow1", hm.IntAttribute(231))
    p1.add_attribute("theta", hm.DoubleAttribute(random.random() * math.pi))
    p1.add_attribute("phi", hm.DoubleAttribute(random.random() * math.pi * 2))

    v2 = hm.GenVertex()
    evt.add_vertex(v2)
    p2 = hm.GenParticle(hm.FourVector(1.0, 1.0, -7000, 7000), 2212, 3)
    evt.add_particle(p2)
    p2.add_attribute("flow1", hm.IntAttribute(243))
    p2.add_attribute("theta", hm.DoubleAttribute(random.random() * math.pi))
    p2.add_attribute("phi", hm.DoubleAttribute(random.random() * math.pi * 2))
    v2.add_particle_in(p2)
    #
    #     create the outgoing particles of v1 and v2
    p3 = hm.GenParticle(hm.FourVector(0.750, -1.569, 32.191, 32.238), 1, 3)
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
    #
    #     create v3
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
    #
    #     create v4
    v4 = hm.GenVertex(hm.FourVector(0.12, -0.3, 0.05, 0.004))
    evt.add_vertex(v4)
    v4.add_particle_in(p5)
    p7 = hm.GenParticle(hm.FourVector(-2.445, 28.816, 6.082, 29.552), 1, 1)
    evt.add_particle(p7)
    v4.add_particle_out(p7)
    p8 = hm.GenParticle(hm.FourVector(3.962, -49.498, -26.687, 56.373), -2, 1)
    evt.add_particle(p8)
    v4.add_particle_out(p8)
    #
    #     tell the event which vertex is the signal process vertex

    evt.add_attribute("signal_process_vertex", hm.IntAttribute(v3.id()))
    print(dir(hm))
    print(hm.Print.content(evt))
    #    we now print it out in old format
    print(hm.Print.listing(evt, 8))
    #     print each particle so we can see the polarization
    for ip in evt.particles():
        print(hm.Print.line(ip, True))

    xout1 = hm.WriterAscii(python_label() + "testBoost1.out")
    xout1.set_precision(6)
    xout1.write_event(evt)
    xout1.close()
    # different outputs
    ff = io.StringIO()
    hm.Print.listing(ff, evt)
    print(ff.getvalue())
    # different outputs
    for ip in evt.particles():
        print(hm.Print.line(ip, True))
    for ip in evt.particles():
        print(hm.Print.line(ip, True))
    xout2 = hm.WriterAscii(python_label() + "testBoost2.out")
    xout2.set_precision(6)
    xout2.write_event(evt)
    xout2.close()
    assert COMPARE_ASCII_FILES(python_label() + "testBoost1.out", python_label() + "testBoost2.out") == 0
    evt.clear()
    return 0


if __name__ == "__main__":
    result = 1
    try:
        result = test_Print()
    except:
        print("FAILED")
        result = 1
    sys.exit(result)
