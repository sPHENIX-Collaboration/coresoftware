# Translated into python using testPolarization.cc (garren@fnal.gov, Oct. 2010) as an example
#
# andrii.verbytskyi@mpp.mpg.org, Nov. 2018
#
from pyHepMC3TestUtils import update_path, python_label, fuse_equal
import sys

sys.path = update_path()

import random, math
from pyHepMC3TestUtils import COMPARE_ASCII_FILES
from pyHepMC3 import HepMC3 as hm
from pyHepMC3 import std as std


def test_Attribute():
    # First create the event container, with Signal Process 20, event number 1
    #
    evt = hm.GenEvent(hm.Units.GEV, hm.Units.MM)
    evt.set_event_number(1)
    evt.add_attribute("signal_process_id", hm.IntAttribute(20))
    x1 = hm.VectorIntAttribute()
    evt.add_attribute("somevector1", x1)
    x2 = hm.VectorDoubleAttribute(std.vector_double([1.0, 2.0, 3.0, 4.0]))
    evt.add_attribute("somevector2", x2)
    v1 = hm.GenVertex()
    evt.add_vertex(v1)
    p1 = hm.GenParticle(hm.FourVector(0, 0, 7000, 7000), 2212, 3)
    evt.add_particle(p1)
    p1.add_attribute("flow1", hm.IntAttribute(231))
    p1.add_attribute("flow1", hm.IntAttribute(233))
    random1 = random.random() * math.pi
    p1.add_attribute("theta", hm.DoubleAttribute(random1))
    random2 = random.random() * math.pi * 2
    p1.add_attribute("phi", hm.DoubleAttribute(random2))

    v2 = hm.GenVertex()
    evt.add_vertex(v2)
    p2 = hm.GenParticle(hm.FourVector(0, 0, -7000, 7000), 2212, 3)
    evt.add_particle(p2)
    p2.add_attribute("flow1", hm.IntAttribute(243))
    random3 = random.random() * math.pi
    p2.add_attribute("theta", hm.DoubleAttribute(random3))
    random4 = random.random() * math.pi * 2
    p2.add_attribute("phi", hm.DoubleAttribute(random4))
    v2.add_particle_in(p2)

    evt_spid_back = hm.IntAttribute()
    evt_spid_back.from_string(evt.attribute("signal_process_id"))
    assert evt_spid_back.value() == 20

    p1_flow1_back = hm.IntAttribute()
    p1_flow1_back.from_string(p1.attribute("flow1"))
    assert p1_flow1_back.value() == 233

    p2_flow1_back = hm.IntAttribute()
    p2_flow1_back.from_string(p2.attribute("flow1"))
    assert p2_flow1_back.value() == 243

    p1_theta_back = hm.DoubleAttribute()
    p1_theta_back.from_string(p1.attribute("theta"))

    assert fuse_equal(p1_theta_back.value(), random1)

    p2_phi_back = hm.DoubleAttribute()
    p2_phi_back.from_string(p2.attribute("phi"))

    assert fuse_equal(p2_phi_back.value(), random4)

    evt.clear()

    return 0


if __name__ == "__main__":
    result = 1
    try:
        result = test_Attribute()
    except:
        result = 1
    sys.exit(result)
