from distutils.core import setup
from Cython.Build   import cythonize
import os
import numpy

os.environ['CFLAGS'] = '-ffast-math'
os.environ['CFLAGS'] = '-O3'        

setup(
      
    ext_modules = cythonize([
                             "./PPICSimulation.pyx", 
                             "./utils.pyx",
                             ],
                            
                            annotate = False,
                            compiler_directives={'language_level' : '3',
                                                 'boundscheck': False,
                                                 'cdivision': True,
                                                 'overflowcheck': False,
                                                 'initializedcheck': False, 
                                                 'wraparound': False
                                                 },
                            ),
    
    include_dirs = [numpy.get_include()]
)
