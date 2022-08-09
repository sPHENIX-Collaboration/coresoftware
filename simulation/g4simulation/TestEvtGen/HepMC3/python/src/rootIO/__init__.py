try:
 import platform
 if platform.system() != 'Darwin':  
   from ctypes import cdll
   libCore = cdll.LoadLibrary("libCore.so")
except:
  print("An exception occurred while loading libCore.so with ctypes")
from .pyHepMC3rootIO import *
