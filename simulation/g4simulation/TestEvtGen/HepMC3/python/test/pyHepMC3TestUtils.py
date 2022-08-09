import re
import sys, os, math


def python_label():
    v = sys.version_info
    impl = ""
    try:
        pypyv = sys.pypy_version_info
        impl = "pypy"
        print(pypyv)
    except:
        pass
    a = str(impl) + str(v[0]) + "." + str(v[1]) + "." + str(v[2])
    print(a)
    return a


def update_path():
    try:
     os.add_dll_directory(os.path.abspath(os.path.join(os.pardir, python_label())))
     os.add_dll_directory(os.getcwd())
     for val in str(os.getenv("PATH")).split(',:;'):
       os.add_dll_directory(os.path.abspath(val))
     for val in str(os.getenv("LD_LIBRARY_PATH")).split(',:;'):
       os.add_dll_directory(os.path.abspath(val))
     for val in str(os.getenv("DYLD_LIBRARY_PATH")).split(',:;'):
       os.add_dll_directory(os.path.abspath(val))
    except:
     pass
    return [os.path.abspath(os.path.join(os.pardir, python_label()))] + [os.getcwd()] + sys.path


def COMPARE_ASCII_FILES(f1, f2):
    file1 = open(f1)
    file2 = open(f2)
    print("Run comparison")
    string1 = " "
    string2 = " "
    j = 0
    while string1 and string2:
        j = j + 1
        string1 = file1.readline()
        string2 = file2.readline()
        if string1 != string2:
            print(j, "-th strings are not equal", "\n")
            print("   ", string1, "\n")
            print("   ", string2, "\n")
            if not re.match(r"HepMC::Version", string1):
                file1.close()
                file2.close()
                return 1
    file1.close()
    file2.close()
    return 0

def fuse_equal(a, b, rtol=0.00001, atol=0.000001):
    if abs(1.0 * a - 1.0 * b) < atol:
        return True
    if abs(1.0 * a - 1.0 * b) < rtol * (abs(1.0 * a + 1.0 * b)):
        return True
    return False
