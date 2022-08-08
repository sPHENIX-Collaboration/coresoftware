from pyHepMC3TestUtils import update_path, python_label
import sys

sys.path = update_path()

from pyHepMC3TestUtils import COMPARE_ASCII_FILES
from pyHepMC3 import HepMC3 as hm


def test_IO1():
    inputA = hm.ReaderAsciiHepMC2("inputIO1.hepmc")
    if inputA.failed():
        sys.exit(1)
    outputA = hm.WriterAscii(python_label() + "frominputIO1.hepmc")
    if outputA.failed():
        sys.exit(2)
    while not inputA.failed():
        evt = hm.GenEvent()
        inputA.read_event(evt)
        if inputA.failed():
            print("End of file reached. Exit.\n")
            break
        outputA.write_event(evt)
        evt.clear()
    inputA.close()
    outputA.close()

    inputB = hm.ReaderAscii(python_label() + "frominputIO1.hepmc")
    if inputB.failed():
        sys.exit(3)
    outputB = hm.WriterAsciiHepMC2(python_label() + "fromfrominputIO1.hepmc")
    if outputB.failed():
        sys.exit(4)
    while not inputB.failed():
        evt = hm.GenEvent()
        inputB.read_event(evt)
        if inputB.failed():
            print("End of file reached. Exit.\n")
            break
        outputB.write_event(evt)
        evt.clear()
    inputB.close()
    outputB.close()
    assert 0 == COMPARE_ASCII_FILES(python_label() + "fromfrominputIO1.hepmc", "inputIO1.hepmc")
    return 0


if __name__ == "__main__":
    result = 1
    try:
        result = test_IO1()
    except:
        result = 1
    sys.exit(result)
