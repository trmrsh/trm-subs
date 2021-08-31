from distutils.core import setup, Extension
from Cython.Build import cythonize
import os, numpy

library_dirs = []
include_dirs = []

# need to direct to where includes and  libraries are
if 'TRM_SOFTWARE' in os.environ:
    library_dirs.append(os.path.join(os.environ['TRM_SOFTWARE'], 'lib64'))
    library_dirs.append(os.path.join(os.environ['TRM_SOFTWARE'], 'lib'))
    include_dirs.append(os.path.join(os.environ['TRM_SOFTWARE'], 'include'))
else:
    print("Environment variable TRM_SOFTWARE pointing to location of shareable libraries and includes not defined!")

include_dirs.append(numpy.get_include())

subs = [Extension("trm.subs._subs",
                 [os.path.join('trm','subs','_subs.pyx')],
                 define_macros = [('MAJOR_VERSION', '1'),
                                  ('MINOR_VERSION', '0')],
                 include_dirs = include_dirs,
                 library_dirs = library_dirs,
                 runtime_library_dirs = library_dirs,
                 libraries = ['subs'])]

setup(name='trm.subs',
      version='1.0.1',
      packages = ['trm', 'trm.subs', 'trm.subs.input', 'trm.subs.plot', 'trm.subs.smtp', 'trm.subs.cpp', 'trm.subs.dvect'],
      ext_modules=cythonize(subs),
      scripts=['scripts/hms2decimal.py'],

      # metadata
      author='Tom Marsh',
      author_email='t.r.marsh@warwick.ac.uk',
      description="Python basic utility module",
      url='http://www.astro.warwick.ac.uk/',
      long_description="""
subs provides an interface to various basic routines as well as a set of routines of general utility.
""",

)

