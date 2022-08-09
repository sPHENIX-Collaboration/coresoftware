from pyHepMC3TestUtils import update_path
import sys

sys.path = update_path()

from pyHepMC3TestUtils import COMPARE_ASCII_FILES
from pyHepMC3 import HepMC3 as hm
import random


def conversion_factor(from_, to_):
    m = hm.FourVector(0.5 - random.random(), 0.5 - random.random(), 0.5 - random.random(), 0.5 - random.random())
    msave = hm.FourVector(m)
    hm.Units.convert(m, from_, to_)
    return m.e() / msave.e()


def neq(a, b):
    if abs(a - b) < 0.001 * (abs(a) + abs(b)):
        return False
    return True


def test_Units():
    err = 0
    evt = hm.GenEvent()

    print("Default hm.Units: ", hm.Units.name(evt.momentum_unit()), " ", hm.Units.name(evt.length_unit()))
    cf = conversion_factor(hm.Units.MomentumUnit.GEV, hm.Units.MomentumUnit.GEV)
    print(cf)
    if neq(cf, 1.0):

        err = err + 1
        print("wrong conversion factor ", cf, " for GEV to GEV - should be 1 \n")

    cf = conversion_factor(hm.Units.MomentumUnit.MEV, hm.Units.MomentumUnit.MEV)
    if neq(cf, 1.0):

        err = err + 1
        print("wrong conversion factor ", cf, " for MEV to MEV - should be 1 \n")

    cf = conversion_factor(hm.Units.MomentumUnit.MEV, hm.Units.MomentumUnit.GEV)
    if neq(cf, 0.001):

        err = err + 1
        print("wrong conversion factor ", cf, " for MEV to GEV - should be 0.001 \n")

    cf = conversion_factor(hm.Units.MomentumUnit.GEV, hm.Units.MomentumUnit.MEV)
    if neq(cf, 1000.0):

        err = err + 1
        print("wrong conversion factor ", cf, " for GEV to MEV - should be 1000 \n")

    cf = conversion_factor(hm.Units.LengthUnit.MM, hm.Units.LengthUnit.MM)
    if neq(cf, 1.0):

        err = err + 1
        print("wrong conversion factor ", cf, " for MM to MM - should be 1 \n")

    cf = conversion_factor(hm.Units.LengthUnit.CM, hm.Units.LengthUnit.CM)
    if neq(cf, 1.0):

        err = err + 1
        print("wrong conversion factor ", cf, " for CM to CM - should be 1 \n")

    cf = conversion_factor(hm.Units.LengthUnit.CM, hm.Units.LengthUnit.MM)
    if neq(cf, 10.0):

        err = err + 1
        print("wrong conversion factor ", cf, " for CM to MM - should be 10 \n")

    cf = conversion_factor(hm.Units.LengthUnit.MM, hm.Units.LengthUnit.CM)
    if neq(cf, 0.1):

        err = err + 1
        print("wrong conversion factor ", cf, " for MM to CM - should be 0.1 \n")

    print(err)
    assert 0 == err
    return 0


if __name__ == "__main__":
    result = 1
    try:
        result = test_Units()
    except:
        result = 1
    sys.exit(result)
