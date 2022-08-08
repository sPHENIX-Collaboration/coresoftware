from pyHepMC3TestUtils import update_path
import sys

sys.path = update_path()

from pyHepMC3TestUtils import COMPARE_ASCII_FILES
from pyHepMC3 import HepMC3 as hm
import random


def test_Pythonization_docs():
    print(hm.GenEvent.particles.__doc__)
    print(hm.GenEvent.event_pos.__doc__)
    print(hm.GenVertex.particles_in.__doc__)
    return 0


if __name__ == "__main__":
    result = 1
    try:
        result = test_Pythonization_docs()
    except:
        result = 1
    sys.exit(result)
