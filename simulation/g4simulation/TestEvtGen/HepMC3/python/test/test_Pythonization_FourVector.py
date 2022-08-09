from pyHepMC3TestUtils import update_path
import sys

sys.path = update_path()

from pyHepMC3TestUtils import COMPARE_ASCII_FILES
from pyHepMC3 import HepMC3 as hm
import random


def test_Pythonization_FourVector():
    evt = hm.GenEvent()
    m = hm.FourVector(0.5 - random.random(), 0.5 - random.random(), 0.5 - random.random(), 0.5 - random.random())
    print(m[0])
    m[1] = 10
    return 0


if __name__ == "__main__":
    result = 1
    try:
        result = test_Pythonization_FourVector()
    except:
        result = 1
    sys.exit(result)
