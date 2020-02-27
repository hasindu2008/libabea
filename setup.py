from distutils.core import setup, Extension
from distutils.sysconfig import *
from distutils.util import *
import os
import os.path
import numpy
from distutils.command.build_py import build_py
from Cython.Distutils import build_ext
from Cython.Build import cythonize


math_lib = ['m']

#extra_compile_args = ['-Wno-strict-prototypes']

py_inc = [get_python_inc()]

np_lib = os.path.dirname(numpy.__file__)
np_inc = [os.path.join(np_lib, 'core/include')]
cmdclass = {'build_py': build_py}

cmdclass.update({'build_ext': build_ext})
packages=['test']

# ext_modules = [Extension("test.abea", 
#                      ["libabea.c", "align.cpp", "events.cpp", "f5c.cpp", "model.cpp",
#                       "abea.pyx"], 
#                       extra_compile_args=["-std=c++11"], 
#                       libraries=math_lib,
#                       include_dirs=py_inc + np_inc)]

extensions = [Extension("test.abea", ["abea.pyx", "libabea.cpp", "align.cpp", "events.cpp", "f5c.cpp", "model.cpp"],
                      extra_compile_args=["-std=c++11"], 
                      libraries=math_lib,
                      include_dirs=py_inc + np_inc)]

setup(name = 'test',
      version='0.0.1',
      requires=['numpy (>=1.3.0)'],
      description='abea',
      author='H',
      author_email='H',
      maintainer='H',
      maintainer_email='H',
      packages=packages,
      cmdclass=cmdclass,
      ext_modules=cythonize(extensions),
      )

