from pyHepMC3TestUtils import update_path
import sys

sys.path = update_path()

from pyHepMC3TestUtils import COMPARE_ASCII_FILES
from pyHepMC3 import HepMC3 as hm
import random


def test_Pythonization_GenEvent():
    evt = hm.GenEvent()
    evt.add_particle(hm.GenParticle())
    print(evt.particles)
    evt.add_vertex(hm.GenVertex())
    return 0


if __name__ == "__main__":
    result = 1
    try:
        result = test_Pythonization_GenEvent()
    except:
        result = 1
    sys.exit(result)
