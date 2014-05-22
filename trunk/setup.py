# -*- coding:utf-8 -*-

from distutils.core import setup
from Cython.Build import cythonize
from setuptools import find_packages


setup(name='BIP',
      version='0.5.14',
      author='Flavio Codeco Coelho',
      author_email='fccoelho@gmail.com',
      url='http://code.google.com/p/bayesian-inference/',
      description='Bayesian Inference Tools for Python',
      zip_safe=False,
      packages=find_packages(),
      #['','BIP','BIP.SDE','BIP.Bayes','BIP.SMC','BIP.Bayes.general','BIP.Bayes.conjugate','BIP.Bayes.tests','BIP.Viz'],
      package_data={'': ['*.txt']},
      # ext_modules=cythonize("BIP/SDE/cgillespie.pyx"),
      install_requires=["numpy", "scipy", "openopt", "liveplots", "cython", "gnuplot-py"],
      test_suite='nose.collector',
      license='GPLv3',
      #        cmdclass = {'build_ext': build_ext},
      #ext_modules=[Extension('BIP/SDE/gillespie', ['BIP/SDE/gillespie.c'])],

)

import os
#This to avoid creating the log file with super-user privileges during instalation.
try:
    os.unlink('/tmp/BIP.log')
except:
    pass
