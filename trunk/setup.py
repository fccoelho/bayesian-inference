# -*- coding:utf-8 -*-

from setuptools import setup
from setuptools import find_packages
from Cython.Build import cythonize
from Cython.Distutils import Extension
from Cython.Distutils import build_ext
import cython_gsl



setup(name='BIP',
      version='0.5.15',
      author='Flavio Codeco Coelho',
      author_email='fccoelho@gmail.com',
      url='http://code.google.com/p/bayesian-inference/',
      description='Bayesian Inference Tools for Python',
      zip_safe=False,
      packages=find_packages(),
      #['','BIP','BIP.SDE','BIP.Bayes','BIP.SMC','BIP.Bayes.general','BIP.Bayes.conjugate','BIP.Bayes.tests','BIP.Viz'],
      package_data={'': ['*.txt']},
      include_dirs=[cython_gsl.get_include()],
      ext_modules=cythonize([Extension("cgillespie", sources=["BIP/SDE/cgillespie.pyx"],
                             libraries=cython_gsl.get_libraries(),
                             library_dirs=[cython_gsl.get_library_dir()],
                             include_dirs=[cython_gsl.get_cython_include_dir()])]),
      install_requires=["numpy", "scipy", "openopt", "liveplots", "cython", "gnuplot-py", "cython_gsl"],
      test_suite='nose.collector',
      license='GPLv3',
      cmdclass={'build_ext': build_ext},
      #ext_modules=[Extension('BIP/SDE/gillespie', ['BIP/SDE/gillespie.c'])],

)

import os
#This to avoid creating the log file with super-user privileges during instalation.
try:
    os.unlink('/tmp/BIP.log')
except:
    pass
