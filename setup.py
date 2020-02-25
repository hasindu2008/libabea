from distutils.core import setup, Extension
from distutils.sysconfig import *
from distutils.util import *
import os
import os.path
import numpy
from distutils.command.build_py import build_py
from Cython.Distutils import build_ext


math_lib = ['m']

extra_compile_args = ['-Wno-strict-prototypes']

py_inc = [get_python_inc()]

np_lib = os.path.dirname(numpy.__file__)
np_inc = [os.path.join(np_lib, 'core/include')]
cmdclass = {'build_py': build_py}

cmdclass.update({'build_ext': build_ext})
packages=['test']

ext_modules = [Extension("test.abea", 
                     ["abea.c", "align.c", "events.c", "f5c.c", "model.c",
                      "abea.pyx"],
                      libraries=math_lib,
                      include_dirs=py_inc + np_inc)]

setup(name = 'test',
      version='0.0.1',
      requires=['numpy (>=1.3.0)'],
      description='dtw2 - dtw algorithm',
      author='H',
      author_email='H',
      maintainer='H',
      maintainer_email='H',
      packages=packages,
      cmdclass=cmdclass,
      ext_modules=ext_modules,
      )

