from pyHepMC3TestUtils import update_path
import sys

sys.path = update_path()

from pyHepMC3TestUtils import COMPARE_ASCII_FILES
from pyHepMC3 import HepMC3 as hm
import random


def test_Pythonization_GenRunInfo():
    ri = hm.GenRunInfo()
    a = hm.GenRunInfo.ToolInfo()
    a.name = "aaa"
    b = hm.GenRunInfo.ToolInfo()
    c = hm.GenRunInfo.ToolInfo()
    b.name = "foo"
    c = ("foo", "1.0", "bar")
    return 0


if __name__ == "__main__":
    result = 1
    try:
        result = test_Pythonization_GenRunInfo()
    except:
        result = 1
    sys.exit(result)
