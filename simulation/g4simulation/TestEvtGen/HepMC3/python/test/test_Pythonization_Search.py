from pyHepMC3TestUtils import update_path
import sys

sys.path = update_path()

from pyHepMC3TestUtils import COMPARE_ASCII_FILES
from pyHepMC3 import HepMC3 as hm
from pyHepMC3.search import HepMC3 as hmsearch
import random, math


def test_Pythonization_Search():
    print(dir(hmsearch))
    inputA = hm.ReaderAsciiHepMC2("inputPythonization_Search.hepmc")
    if inputA.failed():
        print("No input")
        sys.exit(1)
    while not inputA.failed():
        evt = hm.GenEvent()
        inputA.read_event(evt)
        if inputA.failed():
            print("End of file reached. Exit.\n")
            break
        v = evt.vertices()[2]
        evt.clear()
    inputA.close()
    return 0


if __name__ == "__main__":
    result = 1
    try:
        result = test_Pythonization_Search()
    except:
        result = 1
    sys.exit(result)
